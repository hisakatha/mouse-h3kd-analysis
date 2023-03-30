#!/usr/bin/env python3
import pandas as pd
import numpy as np
import pickle
import gzip
import matplotlib.pyplot as plt
import seaborn as sns

prefix = "get_enhancer_ratio_intrachr_1k"
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

# These files must have headers and Geneid column at the first column
clusters = pd.concat([
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_up.tsv", usecols=[0]).assign(cluster="Up"),
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_nochange.tsv", usecols=[0]).assign(cluster="NoChange"),
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_down.tsv", usecols=[0]).assign(cluster="Down"),
]).set_index("Geneid")

samples_1001bp = [
    "2C 17hpi 1001bp embryonic genes",
    "2C control 1001bp embryonic genes",
    "2C H3.1/3.2KD 1001bp embryonic genes",
    ]
#index_1001bp = [samples.index(x) for x in samples_1001bp]
samples_2001bp = [
    "2C 17hpi 2001bp embryonic genes",
    "2C control 2001bp embryonic genes",
    "2C H3.1/3.2KD 2001bp embryonic genes",
    ]
#index_2001bp = [samples.index(x) for x in samples_2001bp]

merged = pd.concat([df.join(clusters).assign(sample=sample) for (sample, df) in feature1_summaries.items()])
merged_1001bp = merged[merged["sample"].isin(samples_1001bp)]
merged_2001bp = merged[merged["sample"].isin(samples_2001bp)]

plt.figure()
sns.boxplot(data=merged_1001bp, x="cluster", y="ratio", hue="sample")
plt.ylabel("Ratio of promoter-enhancer reads to promoter reads")
plt.savefig(prefix + ".boxplot_ratio_by_cluster_1001.svg")

plt.figure()
sns.boxplot(data=merged_2001bp, x="cluster", y="ratio", hue="sample")
plt.ylabel("Ratio of promoter-enhancer reads to promoter reads")
plt.savefig(prefix + ".boxplot_ratio_by_cluster_2001.svg")

plt.figure()
sns.boxplot(data=merged_1001bp, x="cluster", y="num_reads_with_feature2", hue="sample", sym="")
plt.ylabel("#promoter-enhancer reads (w/o outliers)")
plt.savefig(prefix + ".boxplot_num_by_cluster_1001.svg")

plt.figure()
sns.boxplot(data=merged_2001bp, x="cluster", y="num_reads_with_feature2", hue="sample", sym="")
plt.ylabel("#promoter-enhancer reads (w/o outliers)")
plt.savefig(prefix + ".boxplot_num_by_cluster_2001.svg")
