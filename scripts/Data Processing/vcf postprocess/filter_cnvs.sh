#!/bin/bash

for f in *; do
	cd $f
	echo ${f}
	CMD="zcat somatic.vcf.gz | grep -vw '"
	FCT="/chrs_with_cnv/$f"
	CHRS=""
	while read chr; do
		CHRS=$CHRS"$chr\|"
	done <$FCT
	CHRS=${CHRS::-2}
	CMD="$CMD$CHRS'"
	echo "$CMD | wc -l"
	eval $CMD > somatic.nocnv.vcf
	cd ../
done
