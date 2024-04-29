#!/bin/bash

seismic align --cut-nextseq --no-fastqc -o consensus -x fq/virus/tgev-ut_R1.fq.gz -x fq/virus/tgev-ut_R2.fq.gz ref/tgev_refseq.fa

samtools sort consensus/tgev-ut/align/NC_038861.1.bam | samtools consensus -a -f fasta -o ref/tgev_cons.fa -

