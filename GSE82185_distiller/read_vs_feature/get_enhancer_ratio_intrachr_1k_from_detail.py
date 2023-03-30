#!/usr/bin/env python3
# NOTE: This script requires huge RAM
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import gzip
import gc

def load_feature1_summary2(path, doTest=False):
    feature_interaction_detail_dtypes = {
        "pair_name": str,
        "pair_chrom_a": str,
        "pair_start_a": int,
        "pair_end_a": int,
        "pair_score": str,
        "pair_strand_a": str,
        "feature1_chrom": str,
        "feature1_start": int,
        "feature1_end": int,
        "feature1_name": str,
        "feature1_score": str,
        "feature1_strand": str,
        "feature1_overlap_length": int,
        "pair_chrom_b": str,
        "pair_start_b": int,
        "pair_end_b": int,
        "pair_score_duplicate": str,
        "pair_strand_b": str,
        "feature2_chrom": str,
        "feature2_start": int,
        "feature2_end": int,
        "feature2_name": str,
        "feature2_score": str,
        "feature2_strand": str,
        "feature2_overlap_length": int,
    }

    usecols = ["pair_name", "pair_chrom_a", "pair_start_a", "pair_end_a", "feature1_name",
            "pair_chrom_b", "pair_start_b", "pair_end_b", "feature2_name"]
    dd1 = pd.read_table(path, dtype=feature_interaction_detail_dtypes, usecols=usecols, low_memory=True, engine="c")
    print("read_table done. file=" + path)
    dd1["with_feature2"] = (dd1["feature2_name"] != ".")
    dd1["is_intra_chrom"] = (dd1["pair_chrom_a"] == dd1["pair_chrom_b"])
    # abs((start_a + end_a) / 2 - (start_b + end_b) / 2)
    dd1["pair_abs_dist"] = (dd1["pair_start_a"] + dd1["pair_end_a"] - dd1["pair_start_b"] - dd1["pair_end_b"]).abs() / 2

    # Test whether each read pair satisfies conditions
    if doTest:
        dd2 = dd1.groupby(["feature1_name","pair_name"]).agg(num_feature2=("with_feature2", "sum"),
            is_intra_chrom2=("is_intra_chrom", lambda x: len(np.unique(x))),
            pair_abs_dist2=("pair_abs_dist", lambda x: len(np.unique(x))))
        if not all(dd2["is_intra_chrom2"] == 1):
            raise ValueError("Each read pair must have identical chromosome")
        if not all(dd2["pair_abs_dist2"] == 1):
            raise ValueError("Each read pair must have identical distance")

    print("Modification of the initial table (dd1) done. file=" + path)
    dd2 = dd1.groupby(["feature1_name","pair_name"]).agg(num_feature2=("with_feature2", "sum"),
        is_intra_chrom=("is_intra_chrom", "first"),
        pair_abs_dist=("pair_abs_dist", "first"))
    print("groupby (dd2) done. file=" + path)
    # Use read pairs which are intra-chromosomal and distant more than 1kbp
    dd2 = dd2[dd2["is_intra_chrom"] & (dd2["pair_abs_dist"] >= 1000)]
    d_summary = dd2.groupby("feature1_name").agg(num_reads_with_feature2=("num_feature2", lambda x: sum(x > 0)), num_reads=("num_feature2", "size"))
    print("Second groupby (d_summary) done. file=" + path)
    d_summary["ratio"] = d_summary["num_reads_with_feature2"] / d_summary["num_reads"]
    d_summary.drop(".", inplace=True, errors="ignore")
    del dd1
    del dd2
    num_gc = gc.collect()
    print(f"gc.collect() done. return={num_gc}, file={path}")
    return d_summary


paths = {
    "GSE82185 L2C 1001bp PcG targets": "GSE82185_late2C.ext500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.tsv.gz",
    "GSE82185 L2C 1001bp major ZGA": "GSE82185_late2C.ext500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.tsv.gz",
    "GSE82185 L2C 2001bp PcG targets": "GSE82185_late2C.ext1000.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.tsv.gz",
    "GSE82185 L2C 2001bp major ZGA": "GSE82185_late2C.ext1000.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.tsv.gz",
    "GSE82185 L2C 5001bp PcG targets": "GSE82185_late2C.ext2500.bedpe.gz.promoter_PcG_targets_vs_enhancer.feature_interaction_detail.tsv.gz",
    "GSE82185 L2C 5001bp major ZGA": "GSE82185_late2C.ext2500.bedpe.gz.promoter_major_ZGA_vs_enhancer.feature_interaction_detail.tsv.gz",
}

feature1_summaries = {k: load_feature1_summary2(v) for (k,v) in paths.items()}

out_prefix = "get_enhancer_ratio_intrachr_1k"
with gzip.open(out_prefix + ".pkl.gz", "wb") as fp:
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
fig.savefig(out_prefix + ".ratio.svg")

fig, axs = plt.subplots(len(feature1_summaries), sharex=False, sharey=False)
for (i, summary_dict) in enumerate(feature1_summaries.items()):
    (sample, summary) = summary_dict
    ax = axs[i]
    ax.hist(summary["num_reads_with_feature2"], bins=np.linspace(0, 10000, 51), log=True)
    ax.set_title(sample, loc="right", y=0.5)
    #ax.set_ylabel(sample)
fig.supxlabel("#promoter-enhancer reads")
fig.supylabel("Count")
fig.savefig(out_prefix + ".num.svg")

