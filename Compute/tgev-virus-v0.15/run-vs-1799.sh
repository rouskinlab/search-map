#!/bin/bash

set -euxo pipefail


COVS_2KB=../covs-1799
TGEV_2KB=all-noaso
TGEV_2KB_ASO=aso
TGEV_FULL=tgev-full-pool
TGEV_AMPS=tgev_amps_pool

# Align the reads from the TGEV 1.8 kb to the full genome.
if [ ! -f $TGEV_2KB.fq.gz ]
then
	samtools fastq -n -o $TGEV_2KB-long1 $COVS_2KB/out/long1/align/tgev.bam
	samtools fastq -n -o $TGEV_2KB-long2 $COVS_2KB/out/long2/align/tgev.bam
	samtools fastq -n -o $TGEV_2KB-noaso $COVS_2KB/out/noaso/align/tgev.bam
	cat $TGEV_2KB-* > $TGEV_2KB.fq
	rm $TGEV_2KB-*
	gzip $TGEV_2KB.fq
fi
seismic -vv align --no-fastqc --no-cut -y $TGEV_2KB.fq.gz tgev_consensus.fa
if [ ! -f $TGEV_2KB_ASO.fq.gz ]
then
	samtools fastq -n -o $TGEV_2KB_ASO.fq $COVS_2KB/out/aso/align/tgev.bam
	gzip $TGEV_2KB_ASO.fq
fi
seismic -vv align --no-fastqc --no-cut -y $TGEV_2KB_ASO.fq.gz tgev_consensus.fa
seismic -vv relate --batch-size 1200 tgev_consensus.fa out/$TGEV_2KB out/$TGEV_2KB_ASO
seismic -vv mask -s sections-1799.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-finfo-read 0.9 --min-mut-gap 0 out/$TGEV_2KB out/$TGEV_FULL
seismic -vv mask -s sections.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-mut-gap 3 out/$TGEV_2KB out/$TGEV_2KB_ASO
seismic -vv table --no-table-read out/$TGEV_2KB/mask out/$TGEV_2KB_ASO/mask/tgev/lri out/$TGEV_FULL/mask

# Compare the TGEV in-cell data with the 1.8 kb segment.
seismic graph scatter --pdf out/$TGEV_FULL/table/tgev/tgev-1799 out/$TGEV_2KB/table/tgev/tgev-1799
seismic graph scatter --pdf out/$TGEV_AMPS/table/tgev/lri/clust-per-pos.csv out/$TGEV_2KB_ASO/table/tgev/lri/mask-per-pos.csv
seismic graph roc --pdf --struct-file models/fold/tgev/tgev-1799/clip__average.ct out/$TGEV_AMPS/table/tgev/both/clust-per-pos.csv

rm -r log

