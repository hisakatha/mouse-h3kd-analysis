#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip

out_prefix = "get_enhancer_stats_per_promoter_intrachr_1k"
with gzip.open(out_prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

fig, axs = plt.subplots(len(feature1_summaries), sharex=True, sharey=True)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["ratio_max_to_sum"], bins=np.linspace(0, 1, 21), log=False)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Ratio of dominant enhancer's #reads to sum of #reads")
fig.supylabel("Count")
fig.savefig(out_prefix + ".max_read_ratio2.svg")

fig, axs = plt.subplots(len(feature1_summaries), sharex=True, sharey=True)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["num_feature2"], bins=np.linspace(0, 100, 101), log=False)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("#enhancers interacting with each promoter")
fig.supylabel("Count")
fig.savefig(out_prefix + ".num2.svg")

