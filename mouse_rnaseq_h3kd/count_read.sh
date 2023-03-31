#!/usr/bin/env bash
set -ue

TEMPLATE_BED="Aligned.sortedByCoord.out.uniq_collated.template.bed"
cat $TEMPLATE_BED | wc -l > $TEMPLATE_BED.count
cat $TEMPLATE_BED | grep -e "^[1-9XY]\|^MT" | wc -l > $TEMPLATE_BED.count_main
