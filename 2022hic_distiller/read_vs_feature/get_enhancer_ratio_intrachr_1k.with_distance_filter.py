#!/usr/bin/env python3
import sys
import pandas as pd
import pickle
import gzip

def load_feature1_summary3(path, min_dist):
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
        "pair_abs_dist": float,
    }
    table_names = ["feature1_name", "feature2_name", "pair_name", "pair_abs_dist"]
    d_summary = pd.read_table(path, dtype=feature_interaction_detail_dtypes, usecols=table_names, names=table_names)\
            .query("pair_abs_dist >= @min_dist")\
            .assign(with_feature2=lambda df: df.feature2_name != ".")\
            .groupby(["feature1_name","pair_name"]).agg(num_feature2=("with_feature2", "sum"))\
            .groupby("feature1_name").agg(num_reads_with_feature2=("num_feature2", lambda x: sum(x > 0)), num_reads=("num_feature2", "size"))\
            .assign(ratio=lambda df: df.num_reads_with_feature2 / df.num_reads)\
            .drop(".", errors="ignore")
    return d_summary


if len(sys.argv) != 3:
    raise ValueError("#arguments must be 2. Usage: <script> <minimum_distance> <output_prefix>")
min_dist = float(sys.argv[1])
out_prefix = sys.argv[2]

paths = {
    "2C 17hpi 1001bp embryonic genes": "2c_17hpi.ext500.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
    "2C control 1001bp embryonic genes": "2c_control.ext500.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
    "2C H3.1/3.2KD 1001bp embryonic genes": "2c_H3-1_and_3-2KD.ext500.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
    "2C 17hpi 2001bp embryonic genes": "2c_17hpi.ext1000.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
    "2C control 2001bp embryonic genes": "2c_control.ext1000.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
    "2C H3.1/3.2KD 2001bp embryonic genes": "2c_H3-1_and_3-2KD.ext1000.bedpe.gz.promoter_embryonic_vs_enhancer.feature_interaction_detail.intrachr_1k.tsv.gz",
}

feature1_summaries = {k: load_feature1_summary3(v, min_dist) for (k,v) in paths.items()}

with gzip.open(out_prefix + ".pkl.gz", "wb") as fp:
    pickle.dump(feature1_summaries, fp)

