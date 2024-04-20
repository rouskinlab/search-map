#!/bin/bash

set -euxo pipefail


python fse-auc-peaks.py out/tgev-full-pool/graph/tgev/full/aucroll_full__masked_m-ratio-q0_45-9.csv out/tgev-full-pool/graph/tgev/full/aucroll_full__masked_m-ratio-q0_45-9_peaks
python fse-auc-peaks.py $HOME/Desktop/sars2-auc.csv $HOME/Desktop/sars2-auc_peaks

