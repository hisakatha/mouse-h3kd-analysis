#!/usr/bin/env bash
set -ue
export LC_ALL=C
export PATH="/work/00/ga17/share/tools/bedtools2/bin:$PATH"
# Count read pairs overlapping with "feature1",
# and then count the other side of the pairs overlapping with "feature2"

# contact pairs in .bedpe[.gz]
pairs="$1"
# BED features with unique names (with exactly six columns)
feature1="$2"
feature2="$3"
# output prefix
prefix="$4"
sort_opt="-T ."
if [[ ! -e $pairs.region1.bed.gz ]]; then
    zcat -f $pairs | cut -f 1-3,7-8,9 | sort $sort_opt -k 1,1 -k 2,2n | gzip -c > $pairs.region1.bed.gz
fi
if [[ ! -e $pairs.region2.bed.gz ]]; then
    zcat -f $pairs | cut -f 4-6,7-8,10 | sort $sort_opt -k 1,1 -k 2,2n | gzip -c > $pairs.region2.bed.gz
fi
zcat -f $pairs.region1.bed.gz | bedtools intersect -sorted -a - -b $feature1 -wao | gzip -c > $prefix.region1_feature1.bed.gz
zcat -f $pairs.region1.bed.gz | bedtools intersect -sorted -a - -b $feature2 -wao | gzip -c > $prefix.region1_feature2.bed.gz
zcat -f $pairs.region2.bed.gz | bedtools intersect -sorted -a - -b $feature1 -wao | gzip -c > $prefix.region2_feature1.bed.gz
zcat -f $pairs.region2.bed.gz | bedtools intersect -sorted -a - -b $feature2 -wao | gzip -c > $prefix.region2_feature2.bed.gz

# Join by pair_name (4th column)
header="$(echo pair_name,pair_chrom_a,pair_start_a,pair_end_a,pair_score,pair_strand_a,\
feature1_chrom,feature1_start,feature1_end,feature1_name,feature1_score,feature1_strand,feature1_overlap_length,\
pair_chrom_b,pair_start_b,pair_end_b,pair_score_duplicate,pair_strand_b,\
feature2_chrom,feature2_start,feature2_end,feature2_name,feature2_score,feature2_strand,feature2_overlap_length | tr , '\t')"
interactions_detail=$prefix.feature_interaction_detail.tsv.gz
(echo "$header" &&
    join --check-order -t $'\t' -j 4 <(zcat $prefix.region1_feature1.bed.gz | sort $sort_opt -k4,4) \
    <(zcat $prefix.region2_feature2.bed.gz | sort $sort_opt -k4,4) &&
    join --check-order -t $'\t' -j 4 <(zcat $prefix.region2_feature1.bed.gz | sort $sort_opt -k4,4) \
    <(zcat $prefix.region1_feature2.bed.gz | sort $sort_opt -k4,4)) | gzip -c > $interactions_detail

interactions=$prefix.feature_interaction.tsv.gz
# Output: feature1_name, feature2_name, and pair_name delimited by tab
zcat $interactions_detail | tail -n +2 | awk -F "\t" -v OFS="\t" '{print $10,$22,$1}' | sort $sort_opt | uniq | gzip -c > $interactions

