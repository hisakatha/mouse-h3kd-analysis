#!/usr/bin/env bash
set -ue

# Get template (transcribed, mapped/no soft clipped) regions in bed format 
# from a paired-end alignment bam file.
# The TLEN field of an alignment bam file must not include soft clipped regions,
# which is confirmed in the STAR results with "--outSAMtlen 1" (default).
# The alignment bam must have only the primary alignments and be collated by query names
# (using such as "samtools collate -f").
bam="Aligned.sortedByCoord.out.uniq_collated.bam"

# 3rd: RNAME, 4th: POS, 9th: TLEN
# For alignments with positive TLEN, that is, positive strand reads,
# convert 1-based start position into 0-based inclusive start position and exclusive end position
samtools view $bam | awk -v OFS="\t" '$9~/^[1-9][0-9]*$/{print $3,$4-1,$4-1+$9}' > $(basename $bam .bam).template.bed
