#!/usr/bin/env bash
set -ue
input="Aligned.sortedByCoord.out.uniq_collated.template.bed"
# equal to TLEN in sam
cat $input | awk '{print $3 - $2}' > ${input}.length
