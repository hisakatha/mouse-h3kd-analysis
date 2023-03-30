#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
def load_feature1_summary(path):
    d = pd.read_table(path, names=["feature1","feature2","pair"], dtype=str)
    d["withFeature2"] = (d["feature2"] != ".")
    d2 = d.groupby(["feature1","pair"]).agg(numFeature2=("withFeature2", "sum"))
    d_summary = d2.groupby("feature1").agg(numReadsWithFeature2=("numFeature2", lambda x: sum(x > 0)), numReads=("numFeature2", "size"))
    d_summary["ratio"] = d_summary["numReadsWithFeature2"] / d_summary["numReads"]
    d_summary.drop(".", inplace=True)
    return d_summary

paths = {
    "GSE82185 L2C 1001bp PcG targets": "GSE82185_late2C.ext500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction.tsv.gz",
    "GSE82185 L2C 1001bp major ZGA": "GSE82185_late2C.ext500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction.tsv.gz",
    "GSE82185 L2C 2001bp PcG targets": "GSE82185_late2C.ext1000.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction.tsv.gz",
    "GSE82185 L2C 2001bp major ZGA": "GSE82185_late2C.ext1000.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction.tsv.gz",
    "GSE82185 L2C 5001bp PcG targets": "GSE82185_late2C.ext2500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction.tsv.gz",
    "GSE82185 L2C 5001bp major ZGA": "GSE82185_late2C.ext2500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction.tsv.gz",
}

feature1_summaries = {k: load_feature1_summary(v) for (k,v) in paths.items()}

# Save results to a pickle file
import pickle
import gzip
with gzip.open("get_enhancer_ratio.pkl.gz", "wb") as fp:
    pickle.dump(feature1_summaries, fp)

fig, axs = plt.subplots(len(feature1_summaries), sharex=True, sharey=True)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["ratio"], bins=np.linspace(0, 1, 21), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Ratio of promoter-enhancer reads to promoter reads")
fig.supylabel("Count")
fig.savefig("get_enhancer_ratio.ratio.svg")

fig, axs = plt.subplots(len(feature1_summaries), sharex=False, sharey=False)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["numReadsWithFeature2"], bins=np.linspace(0, 10000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("#promoter-enhancer reads")
fig.supylabel("Count")
fig.savefig("get_enhancer_ratio.num.svg")

