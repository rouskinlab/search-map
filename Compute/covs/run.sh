#!/bin/bash

# Compute the Spearman correlation between each pair of +ASO and -ASO profiles
# using a sliding window.
# Predict the structure of each 1799 nt segment based on the mutation rates in
# the -ASO profile.

# For each coronavirus, generate mutational profiles of
# - the 1799 nt (long) segment (2 replicates)
# - the 239 nt (short) segment (2 replicates)
seismic -v align --no-fastqc --min-reads 500 -x fq-size covs.fa
seismic -v relate --min-reads 500 covs.fa out/long* out/short*
seismic -v mask -s sections-fse.csv --min-finfo-read 0.001 --min-ninfo-pos 500 out/long* out/short*
seismic -v table out/long*/mask out/short*/mask

# Compare the mutational profiles for the short and long segments.
short_long_scatter () {
	seismic graph scatter --pdf --no-html out/$2/table/$1/fse/mask-per-pos.csv out/$3/table/$1/fse/mask-per-pos.csv
}
for cov in SARS1 SARS2 BM48 MERS OC43 MHV TGEV B1A IBV CMC;
do
	echo "Short vs. Long: $cov"
	ref=${cov}_1799
	short_long_scatter $ref short1 short2
	short_long_scatter $ref long1 long2
	short_long_scatter $ref short1 long1
	short_long_scatter $ref short2 long2
done

# For each coronavirus, generate mutational profiles of
# - the 1799 nt segment -ASOs (1 replicate)
# - the 1799 nt segment +ASOs (1 replicate)
#seismic -v align --no-fastqc -x fq-aso covs.fa
seismic -v relate covs.fa out/*aso
seismic -v mask --min-finfo-read 0.001 out/*aso
seismic -v table out/*aso/mask
seismic -v fold -q 0.95 covs.fa out/noaso

# Compare the 

