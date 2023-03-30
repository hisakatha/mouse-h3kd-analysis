#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip

# Every pair must be valid.
# That is, every record must not have NA or NULL pair
def load_pair_abs_dist(path, doTest=False):
    table_names = ["pair_name", "pair_abs_dist"]
    table_dtypes = {"pair_name":str, "pair_abs_dist":float}
    d1 = pd.read_table(path, names=table_names, dtype=table_dtypes, usecols=["pair_abs_dist"])
    return d1

paths = {
    "GSE82185 L2C 1001bp PcG targets": "GSE82185_late2C.ext500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
    "GSE82185 L2C 1001bp major ZGA": "GSE82185_late2C.ext500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
    "GSE82185 L2C 2001bp PcG targets": "GSE82185_late2C.ext1000.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
    "GSE82185 L2C 2001bp major ZGA": "GSE82185_late2C.ext1000.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
    "GSE82185 L2C 5001bp PcG targets": "GSE82185_late2C.ext2500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
    "GSE82185 L2C 5001bp major ZGA": "GSE82185_late2C.ext2500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.intrachr_1k.uniq_pairs.tsv.gz",
}

pair_dists = {k: load_pair_abs_dist(v) for (k,v) in paths.items()}

out_prefix = "get_enhancer_promoter_read_stats_intrachr_1k"

fig, axs = plt.subplots(len(pair_dists), sharex=True, sharey=True)
for (i, dist_dict) in enumerate(pair_dists.items()):
    (sample, dist) = dist_dict
    ax = axs[i]
    ax.hist(dist["pair_abs_dist"], bins=np.linspace(0, 200000000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Distance between 5' ends of enhancer-promoter reads")
fig.supylabel("Count")
fig.savefig(out_prefix + ".full.svg")

fig, axs = plt.subplots(len(pair_dists), sharex=True, sharey=True)
for (i, dist_dict) in enumerate(pair_dists.items()):
    (sample, dist) = dist_dict
    ax = axs[i]
    ax.hist(dist["pair_abs_dist"], bins=np.linspace(0, 2000000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Distance between 5' ends of enhancer-promoter reads")
fig.supylabel("Count")
fig.savefig(out_prefix + ".zoom.svg")

fig, axs = plt.subplots(len(pair_dists), sharex=True, sharey=True)
for (i, dist_dict) in enumerate(pair_dists.items()):
    (sample, dist) = dist_dict
    ax = axs[i]
    ax.hist(dist["pair_abs_dist"], bins=np.linspace(0, 100000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("Distance between 5' ends of enhancer-promoter reads")
fig.supylabel("Count")
fig.savefig(out_prefix + ".zoom2.svg")

