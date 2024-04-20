#!/bin/bash

set -euxo pipefail


SAMPLE=sars2_30kb_PCRframeshift

# Process the samples from Lan et al. (2022) over the 2924 nt segment of the SARS-CoV-2 genome.
seismic -vv align --serial -x fq-lan sars2-2924.fa
seismic -vv relate --serial --batch-size 256 sars2-2924.fa out/$SAMPLE[12]
seismic pool -P $SAMPLE out/$SAMPLE[12]
seismic -v mask -s sections-fse.csv --min-finfo-read 0.95 --min-mut-gap 3 out/$SAMPLE*
seismic -vv cluster --serial -k 2 out/$SAMPLE*
seismic -v table out/$SAMPLE*/cluster
# Graph the mutational profiles.
seismic graph profile out/$SAMPLE*/table/sars2-2924/fse/clust-per-pos.csv
# Compare the replicates.
seismic graph scatter out/$SAMPLE[12]/table/sars2-2924/fse/clust-per-pos.csv
# Compare the pooled replicates to the no-ASO control in vitro.
seismic graph scatter --pdf --no-html out/$SAMPLE/table/sars2-2924/fse/clust-per-pos.csv out/ctrl2-deep/table/sars2-2924/fse/clust-per-pos.csv
seismic graph corroll -m pcc -w 25 -n 3 --pdf --no-html out/$SAMPLE/table/sars2-2924/fse/clust-per-pos.csv out/ctrl2-deep/table/sars2-2924/fse/clust-per-pos.csv

rm -r log

