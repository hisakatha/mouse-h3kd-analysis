#!/usr/bin/env python3
import numpy as np
import pickle
import gzip
import matplotlib.pyplot as plt
from scipy import stats

prefix = "get_enhancer_ratio_intrachr_1k_replicate"
with gzip.open(prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

fig, axs = plt.subplots(len(feature1_summaries), sharex=True, sharey=True)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["ratio"], bins=np.linspace(0, 1, 21), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Ratio of promoter-enhancer reads to promoter reads")
fig.supylabel("Count")
fig.savefig(prefix + ".ratio.svg")

fig, axs = plt.subplots(len(feature1_summaries), sharex=False, sharey=False)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["num_reads_with_feature2"], bins=np.linspace(0, 10000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("#promoter-enhancer reads")
fig.supylabel("Count")
fig.savefig(prefix + ".num.svg")

samples = list(feature1_summaries.keys())
ratios = [v["ratio"] for v in feature1_summaries.values()]

plt.figure()
plt.boxplot(ratios, labels=samples)
plt.ylabel("Ratio of promoter-enhancer reads to promoter reads")
plt.xticks(rotation=60)
plt.savefig(prefix + ".boxplot_ratio.svg", bbox_inches='tight', pad_inches=0.1)

nums = [v["num_reads_with_feature2"] for v in feature1_summaries.values()]
plt.figure()
plt.boxplot(nums, labels=samples, sym="")
plt.ylabel("#promoter-enhancer reads (w/o outliers)")
plt.xticks(rotation=60)
plt.savefig(prefix + ".boxplot_num.svg", bbox_inches='tight', pad_inches=0.1)

# Correlation between replicates
nrow = 2
ncol = 3
fig, axs = plt.subplots(nrow, ncol, squeeze=False, figsize=(9,6))
for i in range(nrow):
    for j in range(ncol):
        ax = axs[i, j]
        sample_idx = 2 * (i * ncol + j)
        sample1_key = samples[sample_idx]
        summary1 = feature1_summaries[sample1_key]
        sample2_key = samples[sample_idx + 1]
        summary2 = feature1_summaries[sample2_key]
        joined = summary1.join(summary2, how="outer", lsuffix="_rep1", rsuffix="_rep2")
        #print("1st = " + sample1_key)
        #print(summary1)
        #print("2nd = " + sample2_key)
        #print(summary2)
        #print("Joined:")
        #print(joined)
        #print("Joined with nan:")
        #print(joined.loc[np.isnan(joined.num_reads_rep1)])
        cor = stats.spearmanr(joined["num_reads_with_feature2_rep1"], joined["num_reads_with_feature2_rep2"], nan_policy="omit")
        print(cor)
        joined = joined.loc[(joined["num_reads_with_feature2_rep1"] > 0) & (joined["num_reads_with_feature2_rep2"] > 0)]
        hb = ax.hexbin(joined["num_reads_with_feature2_rep1"], joined["num_reads_with_feature2_rep2"], xscale="log", yscale="log", gridsize=25, bins="log")
        fig.colorbar(hb, ax=ax)
        ax.annotate(f"$\\rho=${cor.statistic:.3f} (p={cor.pvalue:.1e})", xy=(0.05, 0.9), xycoords="axes fraction",
            fontsize="x-small", bbox={"facecolor":"white", "alpha":0.7, "pad":1, "edgecolor":"None"})
        ax.set_xlabel(sample1_key, fontsize="small")
        ax.set_ylabel(sample2_key, fontsize="small")
fig.tight_layout()
fig.savefig(prefix + ".cor_num_interact.png", dpi=300)

fig, axs = plt.subplots(nrow, ncol, squeeze=False, figsize=(9,6))
for i in range(nrow):
    for j in range(ncol):
        ax = axs[i, j]
        sample_idx = 2 * (i * ncol + j)
        sample1_key = samples[sample_idx]
        summary1 = feature1_summaries[sample1_key]
        sample2_key = samples[sample_idx + 1]
        summary2 = feature1_summaries[sample2_key]
        joined = summary1.join(summary2, how="outer", lsuffix="_rep1", rsuffix="_rep2")
        #print("1st = " + sample1_key)
        #print(summary1)
        #print("2nd = " + sample2_key)
        #print(summary2)
        #print("Joined:")
        #print(joined)
        #print("Joined with nan:")
        #print(joined.loc[np.isnan(joined.num_reads_rep1)])
        cor = stats.spearmanr(joined["ratio_rep1"], joined["ratio_rep2"], nan_policy="omit")
        print(cor)
        #joined = joined.loc[(joined["ratio_rep1"] > 0) & (joined["ratio_rep2"] > 0)]
        hb = ax.hexbin(joined["ratio_rep1"], joined["ratio_rep2"], xscale="linear", yscale="linear", gridsize=25, bins="log")
        fig.colorbar(hb, ax=ax)
        ax.annotate(f"$\\rho=${cor.statistic:.3f} (p={cor.pvalue:.1e})", xy=(0.05, 0.9), xycoords="axes fraction",
            fontsize="x-small", bbox={"facecolor":"white", "alpha":0.7, "pad":1, "edgecolor":"None"})
        ax.set_xlabel(sample1_key, fontsize="small")
        ax.set_ylabel(sample2_key, fontsize="small")
fig.tight_layout()
fig.savefig(prefix + ".cor_ratio_interact.png", dpi=300)
