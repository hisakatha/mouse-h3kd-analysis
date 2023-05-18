#!/usr/bin/env bash
set -u
make_cool () {
    local input=$1
    # resolution (bp) for a new file
    local new_res=10000
    local output="$(basename $f).${new_res}.cool"
    if [[ -e $output ]]; then
        echo "Skip make a cool file: $output"
    else
        cooler cp ${input}::resolutions/${new_res} $output
    fi
}

for f in ../results/coolers_library/*.no_filter.1000.mcool; do
    echo $f
    make_cool $f
done
for f in ../../GSE82185_distiller/results/coolers_library/*.no_filter.1000.mcool; do
    echo $f
    make_cool $f
done
for f in ../../GSE80006_distiller_v2/results/coolers_library/*.no_filter.1000.mcool; do
    echo $f
    make_cool $f
done
for f in ../../GSE80006_distiller_v2/results/coolers_library_group/*.no_filter.1000.mcool; do
    echo $f
    make_cool $f
done
