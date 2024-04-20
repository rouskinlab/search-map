#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-00:15
#SBATCH -p short
#SBATCH --array=1-62
#SBATCH -o job_%A_%a.out
#SBATCH -e job_%A_%a.err

python fold_long_fse.py $SLURM_ARRAY_TASK_ID