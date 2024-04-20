#!/bin/bash

set -euxo pipefail

MODELS=models

# Process both replicates against the SARS-2 1.8 kb segment.
seismic -vv align --serial -x fq-lna-1 -x fq-lna-2 sars2-1.8kb.fa
seismic -vv relate --serial --batch-size 256 sars2-1.8kb.fa out
# Pool the replicates.
for lna in $(seq 0 9)
do
	seismic pool -P LNA$lna-pool out/LNA$lna-rep[12]
done
# Mask the individual replicates and their pools.
seismic -v mask -s sections.csv --min-finfo-read 0.95 --min-mut-gap 3 out
# Cluster the -ASO samples and every pool.
seismic -vv cluster --serial -k 8 out/LNA[046789]-rep[12]
seismic -vv cluster --serial -k 6 out/LNA[046789]-pool
seismic -vv cluster --serial -k 2 out/LNA[1235]-pool
seismic -vv table --no-table-read out/*/*s* --serial

# Graph the mutational profiles of the clustered, pooled samples.
seismic graph profile out/LNA?-pool/table/sars2/fse/clust-per-pos.csv
seismic graph profile out/LNA0-rep[12]/table/sars2/fse/clust-per-pos.csv
# Compare every pair of replicates to check overall reproducibility.
for lna in $(seq 0 9)
do
	seismic graph scatter out/LNA$lna-rep[12]/table/sars2/fse/mask-per-pos.csv
done
# Compare the replicates that were clustered to check clustering reproducibility.
for lna in 0 4 6 7 8 9
do
	seismic graph scatter out/LNA$lna-rep[12]/table/sars2/fse/clust-per-pos.csv
done

# Renumber the CT file of the FSE pseudoknot with respect to the 1799 nt segment.
FSE_PK=$MODELS/sars2/pk/zhang2021
if [ ! -f $FSE_PK.ct ]
then
	dot2ct $FSE_PK.db $FSE_PK.ct
	seismic +renumct --inplace -c $FSE_PK.ct 303
fi
# Compute the ROC curve with respect to the pseudoknot for each cluster of samples 0, 4, 6-9.
seismic graph roc --struct-file $FSE_PK.ct out/LNA[046789]-pool/table/sars2/fse/clust-per-pos.csv

# Renumber the CT file of the 1799 nt segment folded with the SARS-2 3kb in vitro data.
FSE_LRI_FILE=fse9__cluster-2-1.ct
FSE_LRI_DIR=$MODELS/sars2/full
seismic +renumct -c ../sars2-2.9kb/out/ctrl2-deep/fold/sars2-2924/sars2-1799/$FSE_LRI_FILE 1 -o $FSE_LRI_DIR
# Compare samples 1, 2, 3, and 5 to the no-ASO condition (LNA0).
for lna in 1 2 3 5
do
	seismic graph corroll -m pcc -w 21 -n 3 --pdf --no-html out/LNA0-pool/table/sars2/fse/mask-per-pos.csv out/LNA$lna-pool/table/sars2/fse/mask-per-pos.csv
	python draw-structure.py $FSE_LRI_DIR/$FSE_LRI_FILE lnas.csv LNA$lna out/LNA0-pool__and__LNA$lna-pool/graph/sars2/fse/corroll_masked_m-ratio-q0_21-3_pcc.csv ../../MainFigures/lnas/lna$lna-correl.pdf
done

rm -r log

