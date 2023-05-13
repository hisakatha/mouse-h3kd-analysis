#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import cooler
import cooltools
import bioframe
from scipy import stats
from cooltools.sandbox import expected_smoothing

def get_chrom_arms(genome_id: str):
    chrom_sizes = bioframe.fetch_chromsizes(genome_id)
    centromeres = bioframe.fetch_centromeres(genome_id)
    chrom_arms = bioframe.make_chromarms(chrom_sizes, centromeres)
    return chrom_arms

# mm10 does not have short arms, so you can omit this call
#chrom_arms = get_chrom_arms("mm10")
chrom_arms = None

resolution = 1000
target_col = "count.avg.smoothed.agg.norm"

def expected_cis_from_path(mcool_path, resolution, chrom_arms):
    target_col_raw = "count.avg.smoothed.agg"
    clr = cooler.Cooler(mcool_path + "::resolutions/" + str(resolution))
    cvd = cooltools.expected_cis(clr, view_df=chrom_arms, nproc=1)
    cvd_cols = {"n_contacts": "count.sum", "contact_freq": "count.avg"}
    cvd_smoothed = expected_smoothing.agg_smooth_cvd(cvd, cols=cvd_cols)
    cvd_smoothed["s_bp"] = cvd_smoothed["dist"] * resolution
    cvd_smoothed_agg = cvd_smoothed.groupby("s_bp", as_index=False).agg({"count.avg.smoothed":"mean"})
    cvd_smoothed_agg.rename(columns={"count.avg.smoothed": target_col_raw}, inplace=True, errors="raise")
    cvd_smoothed_agg.loc[cvd_smoothed_agg["s_bp"] < 2 * resolution, target_col_raw] = np.nan
    base_s_bp = 10000
    base_prob_list = cvd_smoothed_agg.loc[np.isclose(cvd_smoothed_agg["s_bp"], base_s_bp), target_col_raw].to_list()
    if len(base_prob_list) != 1:
        ValueError(f"Input mcool does not have contact info at base {base_s_bp} with resolution {resolution}")
    print(f"Base contact prob. at {base_s_bp} bp = {base_prob_list[0]}")
    cvd_smoothed_agg[target_col] = cvd_smoothed_agg[target_col_raw] / base_prob_list[0]
    print(f"Processed table for {mcool_path}")
    print(cvd_smoothed_agg[:15]); print("..."); print(cvd_smoothed_agg[-5:])
    return cvd_smoothed_agg

clr_path_2c_control = "./coolers_library_group/2c_control.mm10.no_filter.1000.mcool"
contact_freq_2c_control = expected_cis_from_path(clr_path_2c_control, resolution, chrom_arms)
clr_path_2c_17hpi = "./coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool"
contact_freq_2c_17hpi = expected_cis_from_path(clr_path_2c_17hpi, resolution, chrom_arms)
clr_path_2c_h3kd = "./coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool"
contact_freq_2c_h3kd = expected_cis_from_path(clr_path_2c_h3kd, resolution, chrom_arms)

#clr_path_oocyte = "./oocyte_SN_3.mm10.no_filter.1000.mcool"
#contact_freq_oocyte = expected_cis_from_path(clr_path_oocyte, resolution, chrom_arms)

fig, ax = plt.subplots()
ax.loglog(contact_freq_2c_17hpi["s_bp"], contact_freq_2c_17hpi[target_col], "-", label="Early 2-cell")
ax.loglog(contact_freq_2c_control["s_bp"], contact_freq_2c_control[target_col], "-", label="Late 2-cell Control")
ax.loglog(contact_freq_2c_h3kd["s_bp"], contact_freq_2c_h3kd[target_col], "-", label="Late 2-cell H3.1/3.2 KD")
ax.set(
    ylabel="Relative contact frequency",
    xlabel="Distance (bp)",
)
ax.set_aspect(1.0)
ax.legend()
fig.tight_layout()
fig.savefig("contact_prob_no_balancing.svg")

def extract_table_with_tad_length(df):
    min_length = 200_000
    max_length = 1_000_000
    length_col = "s_bp"
    return df[(df[length_col] >= min_length) & (df[length_col] <= max_length)]

contact_freq_2c_17hpi_tad = extract_table_with_tad_length(contact_freq_2c_17hpi)
contact_freq_2c_control_tad = extract_table_with_tad_length(contact_freq_2c_control)
contact_freq_2c_h3kd_tad = extract_table_with_tad_length(contact_freq_2c_h3kd)
prob_col = target_col
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
ax.set_ylim(10**-3, 10**-1)
annot_y1b = 0.06
annot_y1t = annot_y1b * 1.1
annot_y2b = 0.075
annot_y2t = annot_y2b * 1.1
ax.plot([1,1,1.99,1.99], [annot_y1b, annot_y1t, annot_y1t, annot_y1b], c="k")
ax.text(1.5, annot_y1t, f"p={median_test_17hpi_vs_control.pvalue:.1e}", ha="center", va="bottom")
ax.plot([2.01,2.01,3,3], [annot_y1b, annot_y1t, annot_y1t, annot_y1b], c="k")
ax.text(2.5, annot_y1t, f"p={median_test_control_vs_h3kd.pvalue:.1e}", ha="center", va="bottom")
ax.plot([1,1,3,3], [annot_y2b, annot_y2t, annot_y2t, annot_y2b], c="k")
ax.text(2, annot_y2t, f"p={median_test_17hpi_vs_h3kd.pvalue:.1e}", ha="center", va="bottom")
fig.savefig("contact_prob_no_balancing.1mb_boxplot.svg")

fig, ax = plt.subplots(figsize=(3.5,3.5))
ax.loglog(contact_freq_2c_17hpi_tad["s_bp"], contact_freq_2c_17hpi_tad[target_col], "-", label="Early 2-cell")
ax.loglog(contact_freq_2c_control_tad["s_bp"], contact_freq_2c_control_tad[target_col], "-", label="Late 2-cell Control")
ax.loglog(contact_freq_2c_h3kd_tad["s_bp"], contact_freq_2c_h3kd_tad[target_col], "-", label="Late 2-cell H3.1/3.2 KD")
ax.set_ylim(10**-3, 10**-1)
ax.set_xticks(ticks=[], minor=False)
ax.set_xticks(ticks=[200_000, 500_000, 1_000_000], labels=["$2 \\times 10^5$","$5 \\times 10^5$","$10^6$"], minor=True)
ax.set_xlabel("Distance (bp)")
ax.set_ylabel("Relative contact frequency")
ax.legend()
fig.tight_layout()
fig.savefig("contact_prob_no_balancing.1mb.svg")
