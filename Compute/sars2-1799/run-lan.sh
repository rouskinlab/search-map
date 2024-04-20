#!/bin/bash

set -euxo pipefail


SAMPLE=sars2_30kb_PCRframeshift

# Process the in-cell RT-PCR of the FSE from Lan et al. (2022).
seismic -vv align --serial -x fq-lan sars2-1.8kb.fa
seismic -vv relate --serial sars2-1.8kb.fa out/$SAMPLE[12]
seismic pool -P $SAMPLE out/$SAMPLE[12]
seismic -v mask -s sections.csv --min-finfo-read 0.95 --min-mut-gap 3 out/$SAMPLE*
seismic -v cluster -k 6 out/$SAMPLE*
seismic -v table --no-table-read out/$SAMPLE*/cluster
# Compare the clusters of replicates.
seismic graph scatter out/$SAMPLE[12]
# Compare the clusters to the in vitro data.
seismic graph scatter out/$SAMPLE out/LNA0-pool/table/sars2/fse/clust-per-pos.csv

rm -r log

