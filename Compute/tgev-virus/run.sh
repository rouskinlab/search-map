#!/bin/bash

set -euxo pipefail

FULL_POOL=tgev-full-pool
AMPS_POOL=tgev_amps_pool

# Align and relate each sample of total RNA from TGEV-infected cells.
seismic -vv align --serial -x fq-full tgev_consensus.fa
seismic -vv relate --serial --batch-size 1024 tgev_consensus.fa out/tgev-*/align
# Pool the samples of DMS-treated total RNA.
seismic -v pool -P tgev-i-pool out/tgev-i-[12]/relate
seismic -v pool -P tgev-ii-pool out/tgev-ii-[12]/relate
seismic -v pool -P $FULL_POOL out/tgev-i*-pool/relate
# Find the positions in the untreated sample with mutation rates greater than 0.01.
seismic -vv mask -s sections-clip.csv --min-finfo-read 0.9 --min-mut-gap 0 out/tgev-ut
seismic -vv table --no-table-read out/tgev-ut/mask
seismic +listpos --complement --max-fmut-pos 0.01 out/tgev-ut/table/tgev/clip/mask-per-pos.csv
seismic graph profile out/tgev-ut/table/tgev/clip/mask-per-pos.csv
# Mask and tabulate all samples/pools of DMS-treated total RNA.
seismic -vv mask -s sections-clip.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-finfo-read 0.9 --min-mut-gap 0 --min-ninfo-pos 500 out/tgev-i*
seismic -vv mask -s sections-clip.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-finfo-read 0.9 --min-mut-gap 0 out/tgev-full-pool
seismic -vv table --no-table-read out/tgev-i*/mask/tgev/clip
seismic -vv table --no-table-read out/tgev-full-pool/mask/tgev/clip

# Predict the MFE structure of the full TGEV genome with base pairs up to 300 nt.
for fs in 1 1.5 2 2.5 3; do
    seismic -vv fold -s sections-fold$fs.csv -q 0.95 --fold-mfe --fold-md 300 out/$FULL_POOL/table/tgev/clip/mask-per-pos.csv
done
if [ ! -f out/tgev-full-pool/fold/tgev/full/clip__average.db ]
then
    python assemble-tgev-ss.py
fi
# Graph the AUC-ROC between the full genome structure and DMS data.
seismic graph aucroll --struct-sect full --pdf out/$FULL_POOL/table/tgev/clip/mask-per-pos.csv
seismic graph aucroll --struct-file out/$FULL_POOL/fold/tgev/full/clip__average.ct --pdf out/$FULL_POOL/table/tgev/tgev-1799/mask-per-pos.csv
# Graph the mutational profile of the pool.
seismic graph profile out/$FULL_POOL/table/tgev/clip/mask-per-pos.csv
# Compare the mutational profiles of the replicates.
seismic graph scatter --pdf out/tgev-i-[12]/table/tgev/clip/mask-per-pos.csv
seismic graph scatter --pdf out/tgev-ii-[12]/table/tgev/clip/mask-per-pos.csv
seismic graph scatter --pdf out/tgev-i*-pool/table/tgev/clip/mask-per-pos.csv

# For each sample of frameshift stimulating element (FSE) and long-range
# interaction (LRI) amplicons, generate the mutational profile.
seismic -vv align --serial -x fq-amps tgev_consensus.fa
seismic -vv relate --serial --batch-size 4096 tgev_consensus.fa out/tgev_i*
seismic pool -P $AMPS_POOL out/tgev_i*
seismic -v mask -s sections.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-mut-gap 3 out/tgev_*
seismic -vv cluster --serial -k 2 -e 12 out/$AMPS_POOL/mask/tgev/fse out/$AMPS_POOL/mask/tgev/lri
seismic -v join -J both out/tgev_*/mask/tgev/fse out/tgev_*/mask/tgev/lri
seismic -v join -J both -j join_tgev_amps_pool.csv out/$AMPS_POOL/cluster/tgev/fse out/$AMPS_POOL/cluster/tgev/lri
seismic -v table --no-table-read out/tgev_*/mask out/$AMPS_POOL/cluster
seismic graph profile out/$AMPS_POOL/table/tgev/*/*s*-per-pos.csv
seismic graph delprof --compself --no-comppair --pdf out/$AMPS_POOL/table/tgev/*/clust-per-pos.csv
# Predict the structures of the 1,799 nt segment.
seismic -v fold -s sections-1799.csv --fold-percent 100 out/$AMPS_POOL/table/tgev/both/clust-per-pos.csv
seismic graph roc --struct-sect tgev-1799 out/$AMPS_POOL/table/tgev/both/clust-per-pos.csv
seismic graph aucroll --struct-sect tgev-1799 -w 25 -n 5 out/$AMPS_POOL/table/tgev/both/clust-per-pos.csv

# Compare the total RNA and RT-PCR samples.
seismic -v mask -s sections.csv --exclude-file out/tgev-ut/list/tgev/clip/mask-per-pos.csv --min-mut-gap 3 out/tgev-*
seismic -v join -J both out/tgev-*/mask/tgev/fse out/tgev-*/mask/tgev/lri
seismic -v table --no-table-read out/tgev-*/mask/tgev/both
for bio in i ii
do
    for tech in 1 2
    do
	seismic graph scatter --pdf --no-html out/tgev-$bio-$tech/table/tgev/both/mask-per-pos.csv out/tgev_$bio$tech/table/tgev/both/mask-per-pos.csv
    done
done
seismic graph scatter --pdf --no-html out/tgev-full-pool/table/tgev/both/mask-per-pos.csv out/tgev_amps_pool/table/tgev/both/mask-per-pos.csv

rm -r log

