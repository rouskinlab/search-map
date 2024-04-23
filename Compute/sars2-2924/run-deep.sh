#!/bin/bash

set -euxo pipefail


# Process the deep-sequenced samples over the 2924 nt segment of the SARS-CoV-2 genome.
seismic -vv align --serial -x fq-deep sars2-2924.fa
# Extract only the reads aligning to the FSE and amplicon 9.
TEMP=temp
mkdir $TEMP
for sample in ctrl2
do
	BAM=out/$sample-deep/align/sars2-2924.bam
	SAMI=$TEMP/$sample-input.sam
	SAMO=$TEMP/$sample-output.sam
	# Convert the BAM file to a temporary SAM file.
	samtools view -h -o $SAMI $BAM
	# Output only the reads aligning to the FSE/9th amplicon to another SAM file.
	python filter-deep.py $SAMI $SAMO
	# Convert that SAM file back into a BAM file, overwriting the original BAM.
	samtools view -b -o $BAM $SAMO
	rm $SAMI
	rm $SAMO
done
rmdir $TEMP
# Relate the deeply sequenced controls (-ASOs) and +ASOs-9 over the FSE.
seismic -vv relate --serial --batch-size 256 sars2-2924.fa out/*-deep
# Mask everything.
seismic -v mask -s sections-deep.csv --min-finfo-read 0.95 --min-mut-gap 3 out/*-deep
# Cluster control 2.
seismic -vv cluster --serial -k 2 -e 12 --em-thresh 0.01 out/ctrl2-deep/mask/sars2-2924/fse out/ctrl2-deep/mask/sars2-2924/amp9
# Join the FSE and Amplicon 9 into one section, fse9.
seismic -v join -J fse9 out/ctrl2-deep/mask/sars2-2924/fse out/ctrl2-deep/mask/sars2-2924/amp9
seismic -v join -J fse9 -j fse9.csv out/ctrl2-deep/cluster/sars2-2924/fse out/ctrl2-deep/cluster/sars2-2924/amp9
# Tabulate everything.
seismic -vv table out/*-deep
# Compare the control (-ASOs) replicates.
seismic graph scatter --force --pdf --no-html out/ctrl[12]-deep/table/sars2-2924/*/mask-per-pos.csv
seismic graph profile --force --pdf out/ctrl2-deep/table/sars2-2924/*/clust-per-pos.csv
seismic graph histread --force out/ctrl[12]-deep/table/sars2-2924/*/mask-per-read.csv.gz
# Compare the -ASOs and +ASOs-9 samples over the FSE.
seismic graph scatter --force --pdf out/aso9-deep/table/sars2-2924/fse/mask-per-pos.csv out/ctrl2-deep/table/sars2-2924/fse/clust-per-pos.csv
# Compare the -ASOs and +ASOs-4 samples over Amplicon 9.
seismic graph scatter --force --pdf out/aso4-deep/table/sars2-2924/amp9/mask-per-pos.csv out/ctrl2-deep/table/sars2-2924/amp9/clust-per-pos.csv
# Model the structure of the 1,799 nt RNA using the data from the FSE and Amplicon 9.
seismic -vv fold -s sections-1799.csv -q 0.95 --fold-percent 100 --keep-temp out/ctrl2-deep/table/sars2-2924/fse9/clust-per-pos.csv

rm -r log

