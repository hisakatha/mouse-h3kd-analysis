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
pairs="GSE82185_late2C.ext1000.bedpe.gz"
feature1="promoter.major_ZGA_genes.sorted.bed"
feature2="ATAC_early2C_only_DMSO.extend.ensembl.bed"
out_prefix="$pairs.promoter_major_ZGA_vs_enhancer"
./count_interaction.sh $pairs $feature1 $feature2 $out_prefix
