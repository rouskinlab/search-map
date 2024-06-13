#!/bin/bash

### Download results from O2. ###

SRC=ma629@transfer.rc.hms.harvard.edu:/n/data1/hms/microbiology/rouskin/lab/projects/matty/benchmark

rsync -amv --exclude="*.fq" $SRC/sim .
rsync -amv --include="refseq.brickle" --exclude="*.bam" --exclude="*.brickle" --exclude="*.gz" $SRC/out .

