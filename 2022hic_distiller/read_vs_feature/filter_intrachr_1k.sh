#!/usr/bin/env bash
set -u
export LC_ALL=C
sort_opt="-T ."
# input file should be *.feature_interaction_detail.tsv.gz
detail=$1
dist_thres=1000
out_suffix=intrachr_1k
out1="$(basename $detail .tsv.gz).$out_suffix.tsv.gz"
set -o pipefail
zcat $detail | tail -n +2 |
    awk -F "\t" -v dist_thres=$dist_thres -v OFS="\t" \
    '{if(NF!=25){exit 1}; start1=$3; end1=$4; start2=$15; end2=$16; dist=(start1 + end1 - start2 - end2)/2; abs_dist=(dist>0?dist:-dist); chrom1=$2; chrom2=$14; abs_dist=(chrom1==chrom2?abs_dist:-1); if(abs_dist >= dist_thres){feature1=$10; feature2=$22; pair=$1; print feature1,feature2,pair,abs_dist}}' |
    gzip -c > "$out1"

if [[ $? -ne 0 ]]; then
    rm "$out1"
    exit 1
fi

out2=$(basename $out1 .tsv.gz).uniq_pairs.tsv.gz

# Too disk-consuming?
#out3=$(basename $out1 .tsv.gz).uniq_pairs.tsv.gz
#zcat $out1 | tee >(cut -f 1-3 | gzip -c > $out3) | cut -f 3,4 | sort $sort_opt | uniq | gzip -c > $out2

zcat $out1 | cut -f 3,4 | sort $sort_opt | uniq | gzip -c > $out2
