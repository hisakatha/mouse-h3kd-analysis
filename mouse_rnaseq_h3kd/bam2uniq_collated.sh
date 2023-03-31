#!/usr/bin/env bash
set -ue

samtools=samtools
input="Aligned.sortedByCoord.out.bam"
include_flag=0x3 # paired,proper_pair
exclude_flag=0x100 # secondary
$samtools view -b -h -f $include_flag -F $exclude_flag -q 254 -@ 3 $input | $samtools collate -@ 3 -l 6 -f -o $(basename $input .bam).uniq_collated.bam -
