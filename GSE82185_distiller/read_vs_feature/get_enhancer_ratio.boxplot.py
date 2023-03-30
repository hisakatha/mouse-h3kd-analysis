#!/usr/bin/env python3
import pickle
import gzip
import matplotlib.pyplot as plt

with gzip.open("get_enhancer_ratio.pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

samples = list(feature1_summaries.keys())
ratios = [v["ratio"] for v in feature1_summaries.values()]

plt.figure()
plt.boxplot(ratios, labels=samples)
plt.ylabel("Ratio of promoter-enhancer reads to promoter reads")
plt.xticks(rotation=60)
plt.savefig("get_enhancer_ratio.boxplot_ratio.svg", bbox_inches='tight', pad_inches=0.1)

nums = [v["numReadsWithFeature2"] for v in feature1_summaries.values()]
plt.figure()
plt.boxplot(nums, labels=samples, sym="")
plt.ylabel("#promoter-enhancer reads")
plt.xticks(rotation=60)
plt.savefig("get_enhancer_ratio.boxplot_num.svg", bbox_inches='tight', pad_inches=0.1)
