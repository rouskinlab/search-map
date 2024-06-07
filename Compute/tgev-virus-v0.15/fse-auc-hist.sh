#!/bin/bash

set -euxo pipefail


python fse-auc-hist.py out/tgev-full-pool/graph/tgev/full/aucroll_full__masked_m-ratio-q0_45-9.csv out/tgev-full-pool/graph/tgev/full/aucroll_full__masked_m-ratio-q0_45-9_histogram.pdf

