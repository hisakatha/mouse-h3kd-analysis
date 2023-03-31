#!/usr/bin/env bash
set -ue
# Convert raw regions mapped in paired-end data into intergenic bins

RAW_BAM="Aligned.sortedByCoord.out.uniq_collated.bam"
RAW_MERGED="Aligned.sortedByCoord.out.uniq_collated.sorted.merged.bed"
INTERGENIC_BIN="../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene_and_nc.merged.sorted.subtract_bin.bed"
GENIC_BIN="../../mouse_ref/Ensembl_direct/GRCm38/Annotation/Mus_musculus.GRCm38.102.sorted.gene_and_nc.merged.sorted.intersect_bin.bed"
RAW_INTERGENIC_BIN="Aligned.sortedByCoord.out.uniq_collated.sorted.merged.intergene_rounded.bed"
RAW_INTERGENIC_NOOVERHANG="Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.bed"
RAW_INTERGENIC_NOOVERHANG_BIN="Aligned.sortedByCoord.out.uniq_collated.sorted.intergene_nooverhang.merged.rounded.bed"

TMP_RAW_BED="Aligned.sortedByCoord.out.uniq_collated.bed"
TMP_RAW_SORTED="Aligned.sortedByCoord.out.uniq_collated.sorted.bed"
echo "Starting bamtobed"
bedtools bamtobed -i $RAW_BAM > $TMP_RAW_BED
echo "Starting sort"
sort -k 1,1 -k 2,2n $TMP_RAW_BED > $TMP_RAW_SORTED
echo "Starting merge"
bedtools merge -i $TMP_RAW_SORTED > $RAW_MERGED

bedtools subtract -sorted -A -a $TMP_RAW_SORTED -b $GENIC_BIN | bedtools merge > $RAW_INTERGENIC_NOOVERHANG
bedtools intersect -u -sorted -a $INTERGENIC_BIN -b $RAW_INTERGENIC_NOOVERHANG > $RAW_INTERGENIC_NOOVERHANG_BIN

rm $TMP_RAW_BED
rm $TMP_RAW_SORTED
bedtools intersect -u -sorted -a $INTERGENIC_BIN -b $RAW_MERGED > $RAW_INTERGENIC_BIN
