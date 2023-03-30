#!/usr/bin/env python3
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

tpms = pd.read_table("Embryonicgene_5stage_average_tpm.tsv").set_index("Geneid")

prefix = "get_enhancer_ratio_intrachr_1k"
with gzip.open(prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

hic_name_to_long_name = {
    "2C 17hpi 1001bp embryonic genes": "Early 2-cell 1001bp embryonic genes",
    "2C control 1001bp embryonic genes":  "Late 2-cell Control 1001bp embryonic genes",
    "2C H3.1/3.2KD 1001bp embryonic genes": "Late 2-cell H3.1/3.2 KD 1001bp embryonic genes",
    "2C 17hpi 2001bp embryonic genes": "Early 2-cell 2001bp embryonic genes",
    "2C control 2001bp embryonic genes": "Late 2-cell Control 2001bp embryonic genes",
    "2C H3.1/3.2KD 2001bp embryonic genes": "Late 2-cell H3.1/3.2 KD 2001bp embryonic genes",
}
hic_name_to_tpm_name = {
    "2C 17hpi 1001bp embryonic genes": "2c_18hpi",
    "2C control 1001bp embryonic genes": "2c_Control",
    "2C H3.1/3.2KD 1001bp embryonic genes": "2c_H3_1_3_2KD",
    "2C 17hpi 2001bp embryonic genes": "2c_18hpi",
    "2C control 2001bp embryonic genes": "2c_Control",
    "2C H3.1/3.2KD 2001bp embryonic genes": "2c_H3_1_3_2KD",
}

def plot_interaction_vs_expr(min_interaction):
    fig, axs = plt.subplots(len(feature1_summaries), 1, squeeze=False, sharex=True, sharey=True, figsize=(3.5,11))
    for (idx, (sample, feature1_summary)) in enumerate(feature1_summaries.items()):
        ax = axs[idx, 0]
        d2 = feature1_summary.join(tpms)
        sample_name_in_tpm = hic_name_to_tpm_name[sample]
        #ax.scatter(d2["2c_Control"], d2["num_reads_with_feature2"], alpha=0.1)
        #ax.set_xscale("log")
        #ax.set_yscale("log")
        d2 = d2[(d2[sample_name_in_tpm] > 0) & (d2["num_reads_with_feature2"] > min_interaction)]
        hb = ax.hexbin(d2[sample_name_in_tpm], d2["num_reads_with_feature2"], xscale="log", yscale="linear", gridsize=25, extent=(0,4,min_interaction,min_interaction+50), bins="log")
        #hb = ax.hexbin(d2[sample_name_in_tpm], d2["num_reads_with_feature2"], xscale="log", yscale="log", gridsize=25, bins="log")
        ax.set_xticks([1,10**2,10**4])
        ax.set_title(hic_name_to_long_name[sample], fontsize="small")
        ax.set_box_aspect(1)
        cor = stats.spearmanr(d2[sample_name_in_tpm], d2["num_reads_with_feature2"])
        ax.annotate(f"$\\rho=${cor.statistic:.3f} (p={cor.pvalue:.1e})", xy=(0.05, 0.9), xycoords="axes fraction",
            fontsize="x-small", bbox={"facecolor":"white", "alpha":0.7, "pad":1, "edgecolor":"None"})
        fig.colorbar(hb, ax=ax)
    fig.supxlabel("TPM")
    fig.supylabel("#promoter-enhancer reads")
    fig.savefig(prefix + f".plot_interaction{min_interaction}_vs_expr.png", dpi=300)

plot_interaction_vs_expr(0)
plot_interaction_vs_expr(10)

def plot_interaction_ratio_vs_expr(min_interaction=0):
    fig, axs = plt.subplots(len(feature1_summaries), 1, squeeze=False, sharex=True, sharey=True, figsize=(3.5,11))
    for (idx, (sample, feature1_summary)) in enumerate(feature1_summaries.items()):
        ax = axs[idx, 0]
        d2 = feature1_summary.join(tpms)
        sample_name_in_tpm = hic_name_to_tpm_name[sample]
        #ax.scatter(d2["2c_Control"], d2["num_reads_with_feature2"], alpha=0.1)
        d2 = d2[(d2[sample_name_in_tpm] > 0) & (d2["ratio"] >= min_interaction)]
        hb = ax.hexbin(d2[sample_name_in_tpm], d2["ratio"], xscale="log", yscale="linear", gridsize=25, extent=(0,4,min_interaction,1.0), bins="log")
        ax.set_xticks([1,10**2,10**4])
        ax.set_title(hic_name_to_long_name[sample], fontsize="small")
        ax.set_box_aspect(1)
        cor = stats.spearmanr(d2[sample_name_in_tpm], d2["ratio"])
        ax.annotate(f"$\\rho=${cor.statistic:.3f} (p={cor.pvalue:.1e})", xy=(0.05, 0.9), xycoords="axes fraction",
            fontsize="x-small", bbox={"facecolor":"white", "alpha":0.7, "pad":1, "edgecolor":"None"})
        fig.colorbar(hb, ax=ax)
    fig.supxlabel("TPM")
    fig.supylabel("Ratio of #promoter-enhancer reads to #promoter-reads")
    fig.savefig(prefix + f".plot_interaction_ratio{min_interaction}_vs_expr.png", dpi=300)

plot_interaction_ratio_vs_expr()
