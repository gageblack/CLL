SOURCE=$1.vcf.gz
SAMPLES=

## Quality thresholds ##
QUAL_HIGHPASS=20
SAMPLE_DP_HIGHPASS=50
SOMATIC_AO_LOWPASS=5

EXON_BED=GRCh38_reference_genome/GRCh38.gencode.bed.gz
LCR_BED=GRCh38_reference_genome/LCR-hs38.bed.gz
G1K_ACCESSIBILITY_BED=GRCh38_reference_genome/strict_combined.allChr.mask.bed

all: somatic.vcf.gz dp-reports.pdf clean

.PHONY: all

# Required modules:
# vcftools bedtools snpeff vt tabix(samtools)

dp-report.txt: $(addsuffix .DP.FORMAT, $(SAMPLES))
	echo "SAMPLE	DEPTH" > $@
	for sample in $(SAMPLES); do cat $$sample.DP.FORMAT | cut -f3 | grep -v $$sample | while read line; do echo "$$sample	$$line"; done; done >> $@

dp-report500.pdf: dp-report.txt
	echo 'library("ggplot2");' > plot.r
	echo 'data=read.csv("$<", sep="\t", header=T);' >> plot.r
	echo 'p <- ggplot(data, aes(DEPTH, color=SAMPLE)) + geom_density() + xlim(0, 500);' >> plot.r
	echo 'ggsave("$@", p, width=15, height=4, units="in");' >> plot.r
	R CMD BATCH --no-save plot.r
	rm plot.r

dp-report300.pdf: dp-report.txt
	echo 'library("ggplot2");' > plot.r
	echo 'data=read.csv("$<", sep="\t", header=T);' >> plot.r
	echo 'p <- ggplot(data, aes(DEPTH, color=SAMPLE)) + geom_density() + xlim(0, 300);' >> plot.r
	echo 'ggsave("$@", p, width=15, height=4, units="in");' >> plot.r
	R CMD BATCH --no-save plot.r
	rm plot.r

dp-reports.pdf: dp-report300.pdf dp-report500.pdf
	pdfunite dp-report300.pdf dp-report500.pdf $@

## Only use variants found in LCR, G1K accessible, and Exonic regions.
lcr-masked.vcf.gz: $(SOURCE)
	bedtools intersect -a $< -b $(LCR_BED) -wa -v -header | bgzip > $@

masked.vcf.gz: lcr-masked.vcf.gz
	bedtools intersect -a $< -b $(G1K_ACCESSIBILITY_BED) -wa -u -header | bgzip > $@

exonic.vcf.gz: masked.vcf.gz
	bedtools intersect -a $< -b $(EXON_BED) -wa -u -header | bgzip > $@

## Filters variants based on QUALITY score and DEPTH ##
filtered.vcf.gz: exonic.vcf.gz
	zcat $< | snpsift filter "(QUAL > $(QUAL_HIGHPASS)) & (GEN[ALL].DP > $(SAMPLE_DP_HIGHPASS))" | bgzip > $@

biallelic.vcf.gz: filtered.vcf.gz
	vt view -f "N_ALLELE==2" $< -o $@

multiallelic.vcf.gz: filtered.vcf.gz
	vt view -f "N_ALLELE>2" $< -o $@

#SnpEFF annotates and predicts the effects of genetic variants. Keep only somatic variants.
snpeff.vcf.gz: biallelic.vcf.gz
	bash -c "snpeff -i vcf -o vcf -noStats GRCh38.86 $<" | bgzip > $@
	if [ ! $? ]; then rm $@; fi

snpeff-multiallelic.vcf.gz: multiallelic.vcf.gz
	bash -c "snpeff -i vcf -o vcf -noStats GRCh38.86 $<" | bgzip > $@
	if [ ! $? ]; then rm $@; fi

somatic.vcf.gz: snpeff.vcf.gz
	zcat $< | snpsift filter "(GEN[0].AO<$(SOMATIC_AO_LOWPASS))" | bgzip > $@
	if [ ! $? ]; then rm $@; fi

%.DP.FORMAT: biallelic.vcf.gz
	vcftools --gzvcf $< --indv $(basename $(basename $@)) --extract-FORMAT-info DP --out $(basename $(basename $@)) && sed -i 's/\./0/' $@
	if [ ! $? ]; then rm $@; fi

clean:
	mkdir intermediate_files
	mv lcr-masked.vcf.gz masked.vcf.gz exonic.vcf.gz filtered.vcf.gz multiallelic.vcf.gz snpeff-multiallelic.vcf.gz snpeff.vcf.gz *.DP.FORMAT intermediate_files/
	#rm -f exonic.vcf.gz filtered.exonic.vcf.gz snpeff.vcf.gz somatic.vcf.gz *.DP.FORMAT *.log temp-* *.Rout

.PHONY: clean
