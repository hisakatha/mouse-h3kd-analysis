#!/usr/bin/env bash
set -u
pairs="GSE82185_late2C.ext500.bedpe.gz"
feature1="promoter.major_ZGA_genes.sorted.bed"
feature2="ATAC_early2C_only_DMSO.extend.ensembl.bed"
out_prefix="$pairs.promoter_major_ZGA_vs_enhancer"
./count_interaction.sh $pairs $feature1 $feature2 $out_prefix
