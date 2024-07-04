#!/bin/bash

set -euxo pipefail

SCRIPT=plot_genome.py
FIG_DIR=../../MainFigures/search-map
MOD_23S=models/fold/ecoli-23S/full/d.23.b.E.coli.ct
MOD_16S=models/fold/ecoli-16S/full/d.16.b.E.coli.ct

# Process the data with SEISMIC-RNA.
seismic -vv align --serial -x fq ecoli-rrna.fa
seismic -vv relate ecoli-rrna.fa out
seismic -vv mask --min-mut-gap 0 --keep-discontig out/rRNA-nodms-rep1
seismic -vv table out/rRNA-nodms-rep1/mask
seismic +listpos --max-fmut-pos 0.01 --complement out/rRNA-nodms-rep1/table
for ref in ecoli-16S ecoli-23S
do
	seismic -vv mask --mask-pos-file out/rRNA-nodms-rep1/list/$ref/full/mask-per-pos.csv --min-mut-gap 0 --keep-discontig out/*/relate/$ref
done
seismic -vv table out
seismic graph profile out/*/table
seismic graph profile --use-count -r n out/*/table

# Compute the AUC-ROC versus the structure models.
for rna in 16 23
do
	seismic graph roc --struct-file models/fold/ecoli-${rna}S/full/d.$rna.b.E.coli.ct out/rRNA-noaso-rep1/table/ecoli-${rna}S/full/mask-per-pos.csv
	seismic graph aucroll --struct-file models/fold/ecoli-${rna}S/full/d.$rna.b.E.coli.ct out/rRNA-noaso-rep1/table/ecoli-${rna}S/full/mask-per-pos.csv
done

# Graph the rolling correlations with each ASO.
for aso in 0 1a 1b
do
	seismic graph corroll out/rRNA-aso$aso-rep1/table/*/full/mask-per-pos.csv out/rRNA-noaso-rep1/table/*/full/mask-per-pos.csv
	python $SCRIPT models/fold/ecoli-16S/full/d.16.b.E.coli.ct out/rRNA-aso$aso-rep1__and__rRNA-noaso-rep1/graph/ecoli-16S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv $FIG_DIR/aso-$aso.svg
done

rm -r log

