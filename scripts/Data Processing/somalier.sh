#!/bin/bash

module load parallel

REF=GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
SITES=sites.hg38.vcf.gz
echo "Working on extracting each bam file..."


for bam in *.bam; do
	./somalier extract -d extracted/ --sites $SITES -f $REF $bam
	echo "Finished $bam"
done

echo "Done!"
echo "Calculating relatedness..."

./somalier relate -g groups.txt extracted/*.somalier
