#!/usr/bin/env bash
# Activate a conda environment
module load miniconda
source activate /work/ga17/share/tools/conda_envs/basic/

./get_enhancer_ratio_intrachr_1k.with_distance_filter.py 100000 get_enhancer_ratio_intrachr_100k

# Exit from a conda environment
# 'source deactivate' will show warnings
#conda deactivate
