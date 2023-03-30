#!/usr/bin/env bash
# Activate a conda environment
module load miniconda
source activate /work/ga17/share/tools/conda_envs/basic/

./get_enhancer_ratio_intrachr_1k.with_distance_filter.py 1000000 get_enhancer_ratio_intrachr_1m

# Exit from a conda environment
# 'source deactivate' will show warnings
#conda deactivate
