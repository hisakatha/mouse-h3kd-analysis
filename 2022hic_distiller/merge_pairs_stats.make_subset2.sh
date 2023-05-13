#!/usr/bin/env bash
xsv select Sample,total_before_filtering,total,total_mapped,total_nodups,cis,cis_1kb+ pairs_stats.wide.csv > pairs_stats.wide.subset2.csv
