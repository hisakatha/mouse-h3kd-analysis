#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import cooler
import cooltools
import bioframe
from scipy import stats

def get_chrom_arms(genome_id: str):
    chrom_sizes = bioframe.fetch_chromsizes(genome_id)
    centromeres = bioframe.fetch_centromeres(genome_id)
    chrom_arms = bioframe.make_chromarms(chrom_sizes, centromeres)
    return chrom_arms

# mm10 does not have short arms, so you can omit this call
#chrom_arms = get_chrom_arms("mm10")
chrom_arms = None

resolution = 1000

def expected_cis_from_path(mcool_path, resolution, chrom_arms):
    clr = cooler.Cooler(mcool_path + "::resolutions/" + str(resolution))
    cvd = cooltools.expected_cis(clr, view_df=chrom_arms, nproc=1)
    cvd["s_bp"] = cvd["dist"] * resolution
    cvd.loc[cvd["dist"] < 2, "balanced.avg.smoothed.agg"] = np.nan
    cvd_merged = cvd.drop_duplicates(subset=["dist"])[["s_bp", "balanced.avg.smoothed.agg"]]
    return cvd_merged

clr_path_2c_control = "./coolers_library_group/2c_control.mm10.no_filter.1000.mcool"
contact_freq_2c_control = expected_cis_from_path(clr_path_2c_control, resolution, chrom_arms)
clr_path_2c_17hpi = "./coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool"
contact_freq_2c_17hpi = expected_cis_from_path(clr_path_2c_17hpi, resolution, chrom_arms)
clr_path_2c_h3kd = "./coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool"
contact_freq_2c_h3kd = expected_cis_from_path(clr_path_2c_h3kd, resolution, chrom_arms)

#clr_path_oocyte = "./oocyte_SN_3.mm10.no_filter.1000.mcool"
#contact_freq_oocyte = expected_cis_from_path(clr_path_oocyte, resolution, chrom_arms)

fig, ax = plt.subplots()
ax.loglog(contact_freq_2c_17hpi["s_bp"], contact_freq_2c_17hpi["balanced.avg.smoothed.agg"], "-", label="Early 2-cell")
ax.loglog(contact_freq_2c_control["s_bp"], contact_freq_2c_control["balanced.avg.smoothed.agg"], "-", label="Late 2-cell Control")
ax.loglog(contact_freq_2c_h3kd["s_bp"], contact_freq_2c_h3kd["balanced.avg.smoothed.agg"], "-", label="Late 2-cell H3.1/3.2 KD")
#ax.loglog(contact_freq_oocyte["s_bp"], contact_freq_oocyte["balanced.avg.smoothed.agg"], "-", label="oocyte SN#3")
ax.set(
    ylabel="Contact frequency",
    xlabel="Distance (bp)",
)
ax.set_aspect(1.0)
ax.legend()
fig.tight_layout()
fig.savefig("contact_prob.svg")

def extract_table_with_tad_length(df):
    min_length = 200_000
    max_length = 1_000_000
    length_col = "s_bp"
    return df[(df[length_col] >= min_length) & (df[length_col] <= max_length)]

contact_freq_2c_17hpi_tad = extract_table_with_tad_length(contact_freq_2c_17hpi)
contact_freq_2c_control_tad = extract_table_with_tad_length(contact_freq_2c_control)
contact_freq_2c_h3kd_tad = extract_table_with_tad_length(contact_freq_2c_h3kd)
prob_col = "balanced.avg.smoothed.agg"
median_test_17hpi_vs_control = stats.median_test(contact_freq_2c_17hpi_tad[prob_col], contact_freq_2c_control_tad[prob_col])
print("Mood's median test: 17hpi vs control")
print(median_test_17hpi_vs_control)

median_test_17hpi_vs_h3kd = stats.median_test(contact_freq_2c_17hpi_tad[prob_col], contact_freq_2c_h3kd_tad[prob_col])
print("Mood's median test: 17hpi vs H3KD")
print(median_test_17hpi_vs_h3kd)

median_test_control_vs_h3kd = stats.median_test(contact_freq_2c_control_tad[prob_col], contact_freq_2c_h3kd_tad[prob_col])
print("Mood's median test: control vs H3KD")
print(median_test_control_vs_h3kd)

fig, ax = plt.subplots()
ax.boxplot([contact_freq_2c_17hpi_tad[prob_col],
    contact_freq_2c_control_tad[prob_col],
    contact_freq_2c_h3kd_tad[prob_col]],
    labels=["Early 2-cell","Late 2-cell Control","Late 2-cell H3.1/3.2 KD"])
ax.set_ylabel("Contact frequency per 1kb bin within 200kb-1Mb")
ax.set_yscale("log")
ax.set_ylim(10**-5, 10**-3)
annot_y1b = 0.00035
annot_y1t = annot_y1b * 1.1
annot_y2b = 0.00050
annot_y2t = annot_y2b * 1.1
ax.plot([1,1,1.99,1.99], [annot_y1b, annot_y1t, annot_y1t, annot_y1b], c="k")
ax.text(1.5, annot_y1t, f"p={median_test_17hpi_vs_control.pvalue:.1e}", ha="center", va="bottom")
ax.plot([2.01,2.01,3,3], [annot_y1b, annot_y1t, annot_y1t, annot_y1b], c="k")
ax.text(2.5, annot_y1t, f"p={median_test_control_vs_h3kd.pvalue:.1e}", ha="center", va="bottom")
ax.plot([1,1,3,3], [annot_y2b, annot_y2t, annot_y2t, annot_y2b], c="k")
ax.text(2, annot_y2t, f"p={median_test_17hpi_vs_h3kd.pvalue:.1e}", ha="center", va="bottom")
fig.savefig("contact_prob.1mb_boxplot.svg")

fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.loglog(contact_freq_2c_17hpi_tad["s_bp"], contact_freq_2c_17hpi_tad["balanced.avg.smoothed.agg"], "-", label="Early 2-cell")
ax.loglog(contact_freq_2c_control_tad["s_bp"], contact_freq_2c_control_tad["balanced.avg.smoothed.agg"], "-", label="Late 2-cell Control")
ax.loglog(contact_freq_2c_h3kd_tad["s_bp"], contact_freq_2c_h3kd_tad["balanced.avg.smoothed.agg"], "-", label="Late 2-cell H3.1/3.2 KD")
ax.set_ylim(10**-5, 10**-3)
ax.set_xticks(ticks=[], minor=False)
ax.set_xticks(ticks=[200_000, 500_000, 1_000_000], labels=["$2 \\times 10^5$","$5 \\times 10^5$","$10^6$"], minor=True)
ax.set_xlabel("Distance (bp)")
ax.set_ylabel("Contact frequency")
ax.legend()
fig.tight_layout()
fig.savefig("contact_prob.1mb.svg")
