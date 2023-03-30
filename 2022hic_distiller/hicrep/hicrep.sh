#!/usr/bin/env bash
set -u
sample_dir="../results/coolers_library/"
sample_ext=".mm10.no_filter.1000.mcool"
#samples=(2c_17hpi_No1 2c_17hpi_No2)
samples=(2c_17hpi_No1 2c_17hpi_No2 2c_control_No1 2c_control_No2 2c_H3-1_and_3-2KD_No1 2c_H3-1_and_3-2KD_No2)
n_samples=${#samples[@]}
chroms="1 10 11 12 13 14 15 16 17 18 19 2 3 4 5 6 7 8 9 X Y"
for i in $(seq 0 $((n_samples - 1))); do
    for j in $(seq $i $((n_samples - 1))); do
        sample1=${samples[i]}
        sample2=${samples[j]}
        out="hicrep_out.${sample1}.${sample2}"
        echo "Call hicrep $sample1 vs $sample2"
        hicrep --binSize 100000 --h 3 --dBPMax 5000000 \
            ${sample_dir}${sample1}${sample_ext} ${sample_dir}${sample2}${sample_ext} $out \
            --chrNames $chroms
        { echo "chr,hicrep.$sample1.$sample2";
            paste -d, <(echo $chroms | xargs -n1) <(cat $out | grep -v "#" | sed -e "s/ //g"); } > $out.csv
    done
done

