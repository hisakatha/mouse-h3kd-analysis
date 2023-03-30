#!/usr/bin/env bash
# Activate a conda environment
module load miniconda
source activate /work/ga17/share/tools/conda_envs/basic/

./get_enhancer_stats_per_promoter_intrachr_1k.py

# Exit from a conda environment
# 'source deactivate' will show warnings
#conda deactivate
