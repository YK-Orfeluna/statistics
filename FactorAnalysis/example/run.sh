#!/usr/bin/env bash
set -e
# required arguments
inf="BostonHousing_wo_medv.csv"
outd="output"
n_adj=13
# optional arguments
opt="--inf_header --FA_method principal --FA_rotation varimax --kaiser_guttman"
# Run
python $inf $outd $n_adj $opt
exit 0;