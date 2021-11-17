#!/bin/bash

set -v

module load bwa/0.7.17

REF=GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
$ID=
$SM=
$CPUS=

echo 'Sample: '
echo 'ID: '
echo "Reference being used: $REF"

time bwa mem -R '@RG\tID:$ID\tSM:$SM' -t $CPUS $REF $1.fastq.gz $2.fastq.gz | \
	samblaster -r | \
	samtools view -b - > $3.dedup.bam
