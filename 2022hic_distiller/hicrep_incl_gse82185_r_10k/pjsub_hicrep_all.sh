#!/usr/bin/env bash
set -ue
sample_ext=".mm10.no_filter.1000.mcool"
submit_script() {
    local sample1=$1
    local sample2=$2
    sample1_id=$(basename $sample1 $sample_ext)
    sample2_id=$(basename $sample2 $sample_ext)
    script_path=pjsub_hicrep_all.${sample1_id}.${sample2_id}.sh
    if [[ -e $script_path ]]; then
        echo "Skip a job: $sample1 vs $sample2"
        return 0
    fi
    cat << __EOF__ > $script_path
#!/usr/bin/env bash
#PJM -g ga17
#PJM -L node=1
#PJM -L rscgrp=short
#PJM -L elapse=2:00:00
#PJM -S
#PJM -m e,r
#PJM --restart
/work/00/ga17/share/tools/R-install/bin/Rscript hicrep.single.R $sample1_id $sample2_id
__EOF__
    #pjsub $script_path && sleep 60m
    pjsub $script_path && sleep 10m
}

samples=(
../results/coolers_library/2c_17hpi_No1.mm10.no_filter.1000.mcool
../results/coolers_library/2c_17hpi_No2.mm10.no_filter.1000.mcool
../results/coolers_library/2c_control_No1.mm10.no_filter.1000.mcool
../results/coolers_library/2c_control_No2.mm10.no_filter.1000.mcool
../results/coolers_library/2c_H3-1_and_3-2KD_No1.mm10.no_filter.1000.mcool
../results/coolers_library/2c_H3-1_and_3-2KD_No2.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_early2C_No1.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_early2C_No2.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_early2C_No3.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_late2C_No1.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_late2C_No2.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_late2C_No3.mm10.no_filter.1000.mcool
../../GSE82185_distiller/results/coolers_library/GSE82185_late2C_No4.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_3_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_8_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_13_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_14_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_20_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_27_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_29_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_31_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_35_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_37_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_39_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_46_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_51_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_68_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_82_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_83_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_99_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_101_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_114_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_122_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_125_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_130_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_132_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_144_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_150_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_161_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_163_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_165_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_167_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_170_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_172_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_173_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_174_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_175_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_177_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_181_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_185_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_191_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_194_oocyte_SN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_219_oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_1_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_4_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_5_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_10_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_12_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library/GSE80006_15_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_17_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_25_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_30_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_32_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_40_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_41_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_44_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_47_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_48_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_55_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_59_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_62_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_65_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_134_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_152_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_153_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_155_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_159_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_168_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_180_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_182_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_189_oocyte_NSN.mm10.no_filter.1000.mcool
#../../GSE80006_distiller_v2/results/coolers_library/GSE80006_218_oocyte_NSN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library_group/oocyte_SN.mm10.no_filter.1000.mcool
../../GSE80006_distiller_v2/results/coolers_library_group/oocyte_NSN.mm10.no_filter.1000.mcool
)
n_samples=${#samples[@]}

for i in $(seq 0 $((n_samples - 1))); do
    for j in $(seq $((i + 1)) $((n_samples - 1))); do
        submit_script ${samples[i]} ${samples[j]}
    done
done
