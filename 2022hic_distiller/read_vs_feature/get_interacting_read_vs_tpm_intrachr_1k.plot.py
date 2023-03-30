#!/usr/bin/env python3
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt

tpms = pd.read_table("")

prefix = "get_enhancer_ratio_intrachr_1k"
with gzip.open(prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)


