#!/bin/bash

set -euxo pipefail


# Process the tiled samples over the 2924 nt segment of the SARS-CoV-2 genome.
seismic -vv align -x fq-tile sars2-2924.fa
seismic -vv relate --batch-size 1024 sars2-2924.fa out/*-tile/align
seismic mask --min-finfo-read 0.95 --min-mut-gap 3 -s sections-fse.csv out/ctrl?-tile
for n in $(seq 1 12)
do
	seismic mask --min-finfo-read 0.95 --min-mut-gap 3 -s sections-fse.csv out/aso$n-tile
	seismic mask --min-finfo-read 0.95 --min-mut-gap 3 -s sections-amp$n.csv out/aso$n-tile
done
seismic table out/*-tile/mask
# Compute the mutational profiles of each sample over the FSE and the ASO target site.
seismic graph profile --pdf out/*-tile/table/sars2-2924/*/mask-per-pos.csv
# Calculate the correlation between each ensemble average and control replicate 1.
seismic graph corroll -m pcc -w 45 -n 9 --pdf out/ctrl1-tile/table/sars2-2924/fse out/ctrl2-tile/table/sars2-2924/fse
for n in $(seq 1 12)
do
	seismic graph corroll -m pcc -w 45 -n 9 --pdf out/ctrl1-tile/table/sars2-2924/fse out/aso$n-tile/table/sars2-2924/fse
done

rm -r log

