#!/bin/bash

set -euxo pipefail

SCRIPT=../util/pairs_vs_correl.py
FIG_DIR=../../MainFigures/rrna
MOD_23S=models/fold/ecoli-23S/full/d.23.b.E.coli.ct
MOD_16S=models/fold/ecoli-16S/full/d.16.b.E.coli.ct

seismic -vv align -x fq-mv-240125 ecoli-rrna.fa
seismic -vv relate ecoli-rrna.fa out
seismic -v table --no-table-read out/*/relate
# Mask positions with <5000 reads or mutation rates >0.01 in the untreated control.
seismic -vv mask --min-ninfo-pos 5000 --mask-pos ecoli-23S 8 --mask-pos ecoli-23S 174 --mask-pos ecoli-23S 176 --mask-pos ecoli-23S 183 --mask-pos ecoli-23S 184 --mask-pos ecoli-23S 224 --mask-pos ecoli-23S 225 --mask-pos ecoli-23S 227 --mask-pos ecoli-23S 228 --mask-pos ecoli-23S 229 --mask-pos ecoli-23S 264 --mask-pos ecoli-23S 270 --mask-pos ecoli-23S 316 --mask-pos ecoli-23S 441 --mask-pos ecoli-23S 464 --mask-pos ecoli-23S 472 --mask-pos ecoli-23S 542 --mask-pos ecoli-23S 576 --mask-pos ecoli-23S 645 --mask-pos ecoli-23S 646 --mask-pos ecoli-23S 664 --mask-pos ecoli-23S 693 --mask-pos ecoli-23S 694 --mask-pos ecoli-23S 723 --mask-pos ecoli-23S 724 --mask-pos ecoli-23S 727 --mask-pos ecoli-23S 756 --mask-pos ecoli-23S 772 --mask-pos ecoli-23S 773 --mask-pos ecoli-23S 783 --mask-pos ecoli-23S 846 --mask-pos ecoli-23S 852 --mask-pos ecoli-23S 853 --mask-pos ecoli-23S 854 --mask-pos ecoli-23S 973 --mask-pos ecoli-23S 1021 --mask-pos ecoli-23S 1022 --mask-pos ecoli-23S 1023 --mask-pos ecoli-23S 1033 --mask-pos ecoli-23S 1134 --mask-pos ecoli-23S 1135 --mask-pos ecoli-23S 1171 --mask-pos ecoli-23S 1178 --mask-pos ecoli-23S 1188 --mask-pos ecoli-23S 1189 --mask-pos ecoli-23S 1190 --mask-pos ecoli-23S 1211 --mask-pos ecoli-23S 1220 --mask-pos ecoli-23S 1229 --mask-pos ecoli-23S 1420 --mask-pos ecoli-23S 1421 --mask-pos ecoli-23S 1422 --mask-pos ecoli-23S 1423 --mask-pos ecoli-23S 1425 --mask-pos ecoli-23S 1427 --mask-pos ecoli-23S 1434 --mask-pos ecoli-23S 1483 --mask-pos ecoli-23S 1488 --mask-pos ecoli-23S 1553 --mask-pos ecoli-23S 1641 --mask-pos ecoli-23S 1657 --mask-pos ecoli-23S 1658 --mask-pos ecoli-23S 1670 --mask-pos ecoli-23S 1671 --mask-pos ecoli-23S 1693 --mask-pos ecoli-23S 1694 --mask-pos ecoli-23S 1723 --mask-pos ecoli-23S 1725 --mask-pos ecoli-23S 1726 --mask-pos ecoli-23S 1727 --mask-pos ecoli-23S 1730 --mask-pos ecoli-23S 1733 --mask-pos ecoli-23S 1734 --mask-pos ecoli-23S 1735 --mask-pos ecoli-23S 1759 --mask-pos ecoli-23S 1765 --mask-pos ecoli-23S 1789 --mask-pos ecoli-23S 1811 --mask-pos ecoli-23S 1812 --mask-pos ecoli-23S 1865 --mask-pos ecoli-23S 1971 --mask-pos ecoli-23S 2045 --mask-pos ecoli-23S 2060 --mask-pos ecoli-23S 2062 --mask-pos ecoli-23S 2099 --mask-pos ecoli-23S 2132 --mask-pos ecoli-23S 2203 --mask-pos ecoli-23S 2211 --mask-pos ecoli-23S 2233 --mask-pos ecoli-23S 2248 --mask-pos ecoli-23S 2249 --mask-pos ecoli-23S 2250 --mask-pos ecoli-23S 2268 --mask-pos ecoli-23S 2291 --mask-pos ecoli-23S 2292 --mask-pos ecoli-23S 2305 --mask-pos ecoli-23S 2314 --mask-pos ecoli-23S 2441 --mask-pos ecoli-23S 2442 --mask-pos ecoli-23S 2443 --mask-pos ecoli-23S 2507 --mask-pos ecoli-23S 2525 --mask-pos ecoli-23S 2543 --mask-pos ecoli-23S 2573 --mask-pos ecoli-23S 2579 --mask-pos ecoli-23S 2621 --mask-pos ecoli-23S 2622 --mask-pos ecoli-23S 2637 --mask-pos ecoli-23S 2639 --mask-pos ecoli-23S 2640 --mask-pos ecoli-23S 2641 --mask-pos ecoli-23S 2728 --mask-pos ecoli-23S 2750 --mask-pos ecoli-23S 2794 --mask-pos ecoli-23S 2796 --mask-pos ecoli-23S 2799 --mask-pos ecoli-23S 2802 --mask-pos ecoli-23S 2810 --mask-pos ecoli-23S 2820 --mask-pos ecoli-23S 2821 --mask-pos ecoli-23S 2823 --mask-pos ecoli-23S 2824 --mask-pos ecoli-23S 2826 --mask-pos ecoli-23S 2835 --mask-pos ecoli-23S 2836 --mask-pos ecoli-23S 2837 --mask-pos ecoli-23S 2838 --mask-pos ecoli-23S 2839 --mask-pos ecoli-23S 2840 --mask-pos ecoli-23S 2841 --mask-pos ecoli-23S 2853 --mask-pos ecoli-23S 2889 --mask-pos ecoli-23S 2890 --mask-pos ecoli-23S 2891 --mask-pos ecoli-23S 2892 --mask-pos ecoli-23S 2893 --mask-pos ecoli-23S 2897 --mask-pos ecoli-23S 2898 --mask-pos ecoli-23S 2900 --mask-pos ecoli-23S 2903 out/*/relate/ecoli-23S
seismic table --no-table-read out/*/mask
seismic -v graph profile out/*/table/ecoli-*/full/mask-per-pos.csv
seismic -v graph profile --use-count -r n out/*/table/ecoli-*/full/mask-per-pos.csv
seismic -v graph corroll -m pcc out/*/table/ecoli-*/full/mask-per-pos.csv
python $SCRIPT -s $MOD_23S -d out/MV_6_S1_L001__and__MV_7_S2_L001/graph/ecoli-23S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv -a asos/0.csv -o $FIG_DIR/0
python $SCRIPT -s $MOD_23S -d out/MV_6_S1_L001__and__MV_9_S3_L001/graph/ecoli-23S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv -a asos/1a.csv -o $FIG_DIR/1a
python $SCRIPT -s $MOD_23S -d out/MV_10_S4_L001__and__MV_6_S1_L001/graph/ecoli-23S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv -a asos/1b.csv -o $FIG_DIR/1b
python $SCRIPT -s $MOD_23S -d out/MV_11_S5_L001__and__MV_6_S1_L001/graph/ecoli-23S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv -a asos/2a.csv -o $FIG_DIR/2a
python $SCRIPT -s $MOD_23S -d out/MV_12_S6_L001__and__MV_6_S1_L001/graph/ecoli-23S/full/corroll_masked_m-ratio-q0_45-9_pcc.csv -a asos/2b.csv -o $FIG_DIR/2b

rm -r log

