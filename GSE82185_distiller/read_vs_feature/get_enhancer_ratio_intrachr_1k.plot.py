#!/usr/bin/env python3
import pickle
import gzip
import matplotlib.pyplot as plt

prefix = "get_enhancer_ratio_intrachr_1k"
with gzip.open(prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

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

