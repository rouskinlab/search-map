#!/bin/bash

set -euxo pipefail


OUT=out
DEST=../../MainFigures/covs
PVC=../util/pairs_vs_correl.py


#rm -r $OUT/*/align
#rm -r $OUT/*/relate
#rm -r $OUT/*/mask
#rm -r $OUT/*/table
#rm -r $OUT/*/fold
#rm -r $OUT/*/graph
#rmdir $OUT/*__and__*


# LONG VS SHORT

# - 1799 nt (long) segment (2 replicates)
# - 239 nt (short) segment (2 replicates)
seismic align -N 512 -o $OUT -x fq-size covs.fa
seismic relate -N 512 -o $OUT covs.fa $OUT/long[12] $OUT/short[12]

# Compare replicates:
# - long 1 vs long 2
# - short 1 vs short 2
seismic mask -s sections-fse.csv --min-finfo-read 0.05 --min-ninfo-pos 512 $OUT/long[12] $OUT/short[12]
seismic table $OUT/*[12]/mask/*/fse
seismic graph scatter --pdf --no-html -o $OUT $OUT/long[12]
seismic graph scatter --pdf --no-html -o $OUT $OUT/short[12]

# Combine replicates into pooled samples
seismic pool -P both-long $OUT/long[12]
seismic pool -P both-short $OUT/short[12]
seismic mask -s sections-fse.csv --min-finfo-read 0.05 --min-ninfo-pos 1024 $OUT/both-long $OUT/both-short
seismic table $OUT/both-long/mask/*/fse $OUT/both-short/mask/*/fse

# Contrast the long and short mutational profiles
seismic graph scatter --pdf --no-html -o $OUT $OUT/both-long $OUT/both-short


# LONG -ASO VS +ASO

# - 1799 nt (long) -ASOs (1 replicate)
# - 1799 nt (long) +ASOs (1 replicate)
seismic align -N 0 -o $OUT -x fq-aso covs.fa
seismic -qq relate --serial -N 1024 -o $OUT covs.fa $OUT/aso $OUT/noaso

# Compare the -ASO sample to the pool of both long replicates
seismic mask -s sections-clipped.csv --min-finfo-read 0.01 --min-ninfo-pos 512 $OUT/noaso $OUT/both-long
seismic table $OUT/noaso/mask/*/clip $OUT/both-long/mask/*/clip
seismic graph corroll -m scc -w 45 -n 9 --no-html -o $OUT $OUT/noaso $OUT/both-long

# Combine the -ASO sample and both long replicates
seismic pool -P all-noaso $OUT/noaso $OUT/both-long

# Compare the -ASO and +ASO samples
seismic mask -s sections-clipped.csv --min-finfo-read 0.01 --min-ninfo-pos 512 $OUT/aso $OUT/all-noaso
seismic table $OUT/aso/mask/*/clip $OUT/all-noaso/mask/*/clip
seismic graph profile --pdf --no-html -o $OUT $OUT/aso $OUT/all-noaso
seismic graph corroll -m scc -w 45 -n 9 --pdf --no-html -o $OUT $OUT/aso $OUT/all-noaso

# Also predict the RNA structures and compute ROC curves
seismic fold -q 0.95 --fold-max 20 --fold-percent 100 $OUT/all-noaso
for VIRUS in b1a bm48 cmc ibv mhv sars1 sars2 tgev
do
	seismic graph roc --pdf --no-html -o $OUT --struct-file $OUT/all-noaso/fold/$VIRUS/full/clip__average.ct $OUT/all-noaso/table/$VIRUS/clip
done

OUT_REP=$DEST/rep
OUT_ASO=$DEST/aso
mkdir -p $OUT_REP
mkdir -p $OUT_ASO
for VIRUS in b1a bm48 cmc ibv mhv sars1 sars2 tgev
do
	STRUCT=$OUT/all-noaso/fold/$VIRUS/full/clip__average.ct
	ASO_COORDS=$DEST/aso-coords.csv
	NOASO_REPS_CORR=$OUT/both-long__and__noaso/graph/$VIRUS/clip/corroll_masked_m-ratio-q0_45-9_scc.csv
	ASO_VS_NOASO_CORR=$OUT/all-noaso__and__aso/graph/$VIRUS/clip/corroll_masked_m-ratio-q0_45-9_scc.csv
	$PVC -s $STRUCT -d $NOASO_REPS_CORR -a $ASO_COORDS -o $OUT_REP
	$PVC -s $STRUCT -d $ASO_VS_NOASO_CORR -a $ASO_COORDS -o $OUT_ASO
done

# Delete log files if the above commands are all successful.
rm -r log

