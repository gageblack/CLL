#!/bin/bash

set -v
module load freebayes vcflib htslib parallel

REF=GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
REGIONS=regions.10mb.txt
PATIENT=
SAMPLES=
$CPUS=

time freebayes-parallel \
    $REGIONS \
    $CPUS \
    --min-alternate-fraction 0.05 \
    --allele-balance-priors-off \
    --report-genotype-likelihood-max \
    --genotype-qualities \
    --pooled-discrete \
    --pooled-continuous \
    -f $REF \
    -s $SAMPLES \
    -b $1.processed.bam \
    -b $2.processed.bam \
    -b $3.processed.bam \
    -b $4.processed.bam \
    > $PATIENT.vcf
