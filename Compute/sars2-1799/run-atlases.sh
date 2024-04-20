#!/bin/bash

# Make atlas plots over the pseudoknot section (positions 303-371 of the 1799 nt segment) +/- ASOs.

set -euxo pipefail

SCRIPT="atlas-plot.py"
MODELS=models/sars2/mixed/mixed.ct
OUT="../../MainFigures/lnas"


python $SCRIPT $OUT/atlas-LNA0 $MODELS out/LNA0-pool/table/sars2/mixed/clust-per-pos.csv out/LNA0-pool/table/sars2/fse/clust-freq.csv 6

python $SCRIPT $OUT/atlas-LNA4 $MODELS out/LNA4-pool/table/sars2/mixed/clust-per-pos.csv out/LNA4-pool/table/sars2/fse/clust-freq.csv 4

python $SCRIPT $OUT/atlas-LNA6 $MODELS out/LNA6-pool/table/sars2/mixed/clust-per-pos.csv out/LNA6-pool/table/sars2/fse/clust-freq.csv 6

python $SCRIPT $OUT/atlas-LNA7 $MODELS out/LNA7-pool/table/sars2/mixed/clust-per-pos.csv out/LNA7-pool/table/sars2/fse/clust-freq.csv 4

python $SCRIPT $OUT/atlas-LNA8 $MODELS out/LNA8-pool/table/sars2/mixed/clust-per-pos.csv out/LNA8-pool/table/sars2/fse/clust-freq.csv 3

python $SCRIPT $OUT/atlas-LNA9 $MODELS out/LNA9-pool/table/sars2/mixed/clust-per-pos.csv out/LNA9-pool/table/sars2/fse/clust-freq.csv 4

python $SCRIPT $OUT/atlas-cell $MODELS out/sars2_30kb_PCRframeshift/table/sars2/mixed/clust-per-pos.csv out/sars2_30kb_PCRframeshift/table/sars2/fse/clust-freq.csv 4

