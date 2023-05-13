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
    plotpup.py --input_pups $pileup --plot_ticks --no_score --height 2.5 --scale log --vmax 1.3 -o $pileup.pdf
}

make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.head100k.bedpe test_coolpuppy.orig.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.head100k.shift_r1_50k.bedpe test_coolpuppy.shift_r1.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.head100k.shift_r2_50k.bedpe test_coolpuppy.shift_r2.clpy

make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.shift_r1_50k.bedpe test_coolpuppy.full_shift_r1.clpy
make_pileup_plot coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool ../annotation/enhancer_tss_pairs.embryonic_up.shift_r2_50k.bedpe test_coolpuppy.full_shift_r2.clpy
