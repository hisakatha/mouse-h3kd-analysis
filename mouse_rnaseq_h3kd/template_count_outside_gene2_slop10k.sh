#!/usr/bin/env bash
set -ue
# Count read-pairs outside slopped/extended gene regions

TEMPLATE_BED="Aligned.sortedByCoord.out.uniq_collated.template.bed"
TEMPLATE_SORTED="Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed"
WIDTH="10k"
GENIC_REGION="../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene.merged.sorted.slop$WIDTH.merged.bed"

# output files
TEMPLATE_INTERGENIC_NOOVERHANG_READ="Aligned.sortedByCoord.out.uniq_collated.template.sorted.nooverhang_outside_gene2_slop$WIDTH.bed"

if [[ ! -f $TEMPLATE_SORTED ]]; then
    echo "Starting sort"
    sort -k 1,1 -k 2,2n $TEMPLATE_BED > $TEMPLATE_SORTED
fi

bedtools subtract -sorted -A -a $TEMPLATE_SORTED -b $GENIC_REGION > $TEMPLATE_INTERGENIC_NOOVERHANG_READ
cat $TEMPLATE_INTERGENIC_NOOVERHANG_READ | wc -l > $TEMPLATE_INTERGENIC_NOOVERHANG_READ.count
cat $TEMPLATE_INTERGENIC_NOOVERHANG_READ | grep -c -e "^[1-9XY]\|^MT" > $TEMPLATE_INTERGENIC_NOOVERHANG_READ.count_main || [[ $? -eq 1 ]]
