#!/usr/bin/env bash
relpath=..
$relpath/"bam2uniq_collated.sh"
# For template region bin
$relpath/"bam_proper_pair_to_template_region_bed.sh"
$relpath/"template_region_to_rounded_bin.sh"
# For raw region bin
$relpath/"raw_region_to_rounded_bin.sh"
