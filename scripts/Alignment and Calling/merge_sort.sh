#!/bin/bash

set -v

module load bwa/0.7.17

$SAMPLE=
$CPUS=

echo 'Sample: $SAMPLE'

time samtools merge -u -@ $CPUS - ${SAMPLE}1.dedup.bam ${SAMPLE}2.dedup.bam ${SAMPLE}3.dedup.bam | \
	samtools sort -n -@ $CPUS - | \
	samtools view -h - | \
	samblaster -r | \
	samtools sort -@ $CPUS -o ${SAMPLE}.processed.bam

