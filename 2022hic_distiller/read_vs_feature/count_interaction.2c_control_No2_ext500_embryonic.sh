#!/usr/bin/env bash
#PJM -g ga17
#PJM -L node=1
#PJM -L rscgrp=regular
#PJM -L elapse=12:00:00
# #PJM -L rscgrp=short
#PJM -S
#PJM -m e,r
#PJM --restart
set -u
pairs="2c_control_No2.ext500.bedpe.gz"
feature1="promoter.embryonic_igenomes_ucsc.ensembl.bed"
feature2="ATAC_early2C_only_DMSO.extend.ensembl.bed"
out_prefix="$pairs.promoter_embryonic_vs_enhancer"
./count_interaction.sh $pairs $feature1 $feature2 $out_prefix
