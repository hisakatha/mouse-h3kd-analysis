#!/usr/bin/env bash
make_pileup_plot() {
    # input mcool file
    local mcool=$1
    # input feature pair file in BEDPE format
    local feature_pairs=$2
    # output pile-up name
    local pileup=$3
    if [[ ! -e $pileup ]]; then
        coolpup.py ${mcool}::resolutions/10000 $feature_pairs --features_format bedpe --nshifts 1 -o $pileup
    fi
    #plotpup.py --input_pups $pileup --no_score --height 2.5 --scale linear --not_symmetric --vmin 0.6 --vmax 1.4 -o $pileup.pdf
    plotpup.py --input_pups $pileup --no_score --height 2.5 --scale log --vmax 1.3 -o $pileup.pdf
}

make_pileup_plot coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.bedpe 2c_17hpi.embryonic_up.clpy
make_pileup_plot coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_nochange.bedpe 2c_17hpi.embryonic_nochange.clpy
make_pileup_plot coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_down.bedpe 2c_17hpi.embryonic_down.clpy
make_pileup_plot coolers_library_group/2c_control.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.bedpe 2c_control.embryonic_up.clpy
make_pileup_plot coolers_library_group/2c_control.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_nochange.bedpe 2c_control.embryonic_nochange.clpy
make_pileup_plot coolers_library_group/2c_control.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_down.bedpe 2c_control.embryonic_down.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.bedpe 2c_H3-1_and_3-2KD.embryonic_up.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_nochange.bedpe 2c_H3-1_and_3-2KD.embryonic_nochange.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_down.bedpe 2c_H3-1_and_3-2KD.embryonic_down.clpy
