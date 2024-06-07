#!/bin/bash

set -euxo pipefail

seismic -vv align --serial -x fq ecoli-rrna.fa
seismic -vv relate --serial --batch-size 1000000 ecoli-rrna.fa out/*/align
seismic -vv mask --min-ninfo-pos 500 out/*/relate
seismic -vv table --no-table-read out
seismic graph profile out/*/table
seismic graph profile --use-count -r n out/*/table

rm -r log

