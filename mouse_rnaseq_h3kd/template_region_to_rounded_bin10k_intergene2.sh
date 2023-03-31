#!/usr/bin/env bash
set -ue
# Convert raw regions mapped in paired-end data into intergenic bins

TEMPLATE_BED="Aligned.sortedByCoord.out.uniq_collated.template.bed"
TEMPLATE_SORTED="Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed"
TEMPLATE_MERGED="Aligned.sortedByCoord.out.uniq_collated.template.sorted.merged.bed"
WIDTH="10k"
INTERGENIC_BIN="../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene.merged.sorted.subtract_bin$WIDTH.bed"
GENIC_BIN="../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene.merged.sorted.intersect_bin$WIDTH.bed"
# output files
TEMPLATE_INTERGENIC_BIN="Aligned.sortedByCoord.out.uniq_collated.template.sorted.merged.intergene2_rounded_bin$WIDTH.bed"
TEMPLATE_INTERGENIC_NOOVERHANG_READ="Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin$WIDTH.bed"
TEMPLATE_INTERGENIC_NOOVERHANG="Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang_bin$WIDTH.merged.bed"
TEMPLATE_INTERGENIC_NOOVERHANG_BIN="Aligned.sortedByCoord.out.uniq_collated.template.sorted.intergene2_nooverhang.merged.rounded_bin$WIDTH.bed"

if [[ ! -f $TEMPLATE_SORTED ]]; then
    echo "Starting sort"
    sort -k 1,1 -k 2,2n $TEMPLATE_BED > $TEMPLATE_SORTED
fi
if [[ ! -f $TEMPLATE_MERGED ]]; then
    echo "Starting merge"
    bedtools merge -i $TEMPLATE_SORTED > $TEMPLATE_MERGED
fi

bedtools subtract -sorted -A -a $TEMPLATE_SORTED -b $GENIC_BIN > $TEMPLATE_INTERGENIC_NOOVERHANG_READ
cat $TEMPLATE_INTERGENIC_NOOVERHANG_READ | wc -l > $TEMPLATE_INTERGENIC_NOOVERHANG_READ.count
cat $TEMPLATE_INTERGENIC_NOOVERHANG_READ | grep -e "^[1-9XY]\|^MT" | wc -l > $TEMPLATE_INTERGENIC_NOOVERHANG_READ.count_main

bedtools merge -i $TEMPLATE_INTERGENIC_NOOVERHANG_READ > $TEMPLATE_INTERGENIC_NOOVERHANG
bedtools intersect -c -sorted -a $INTERGENIC_BIN -b $TEMPLATE_INTERGENIC_NOOVERHANG_READ > $TEMPLATE_INTERGENIC_NOOVERHANG_BIN

#rm $TEMPLATE_SORTED
bedtools intersect -c -sorted -a $INTERGENIC_BIN -b $TEMPLATE_SORTED > $TEMPLATE_INTERGENIC_BIN
