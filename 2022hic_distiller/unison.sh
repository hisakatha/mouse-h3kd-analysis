#!/usr/bin/env bash
set -u
#root="obcx:/work/ga17/share/2022hic_distiller/./results"
root="ssh://obcx05//work/ga17/share/2022hic_distiller"
unison $* results $root/results\
    -ignore "Path ?*"\
    -ignorenot "Path fastqc/*.html"\
    -ignorenot "Path mapped_parsed_sorted_chunks/*.fastp.{html,json}"\
    -ignorenot "Path pairs_library/*.stats"\
    -ignorenot "Path stats_library_group/*.stats"

#"$root/coolers_library_group/*.no_filter.*.mcool"

#root="obcx:/work/ga17/share/2022hic_distiller/./read_vs_feature"
#rsync -avR $* "$root/get_enhancer_*" ./
#unison $* read_vs_feature $root/read_vs_feature -path "get_enhancer_*"
unison $* read_vs_feature $root/read_vs_feature -ignore "Name ?*" -ignorenot "Name get_enhancer_*"\
    -ignorenot "Name coolbox_plot.*" -ignorenot "Name promoter.embryonic_igenomes_ucsc.ensembl.bed"\
    -ignorenot "Name select_genes_by_tpm_change.*"
unison $* hicrep $root/hicrep
unison $* hicrep_r $root/hicrep_r -ignore "Name *.cool"
unison $* hicrep_r_10k $root/hicrep_r_10k -ignore "Name *.cool"
unison $* hicrep_incl_gse82185 $root/hicrep_incl_gse82185 -ignore "Name *.cool"
unison $* hicrep_incl_gse82185_r_10k $root/hicrep_incl_gse82185_r_10k -ignore "Name *.cool"
