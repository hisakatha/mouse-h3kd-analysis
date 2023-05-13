#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
from functools import reduce
import json
samples = ["2c_17hpi_No1", "2c_17hpi_No2", "2c_control_No1", "2c_control_No2", "2c_H3-1_and_3-2KD_No1", "2c_H3-1_and_3-2KD_No2",]
files_prefix = "./pairs_library/"
files_suffix = ".mm10.dedup.stats"
files = {k: files_prefix + k + files_suffix for k in samples}
files["oocyte_SN_3"] = "../test_Flyamer2017/SRR3343952/oocyte_SN_3.mm10.dedup.stats"
tables = {k:pd.read_table(v, names=["Item", "Count"]).assign(Sample=k) for (k,v) in files.items()}
merged = pd.concat(tables.values())

merged2 = merged.pivot(index="Item", columns="Sample", values="Count")
merged2.to_csv("pairs_stats.csv")

merged_wide = merged.pivot(index="Sample", columns="Item", values="Count")

def get_num_reads_before_filtering(path):
    with open(path) as f:
        return json.load(f)["summary"]["before_filtering"]["total_reads"] / 2

fastp_prefix = "./mapped_parsed_sorted_chunks/"
fastp_suffix = ".lane1.mm10.0.fastp.json"
fastp_files = {k: fastp_prefix + k + fastp_suffix for k in samples}
fastp_files["oocyte_SN_3"] = "../test_Flyamer2017/SRR3343952/oocyte_SN_3.lane1.mm10.0.fastp.json"
num_reads_before_filtering = {k: get_num_reads_before_filtering(v) for (k,v) in fastp_files.items()}
merged_wide["total_before_filtering"] = num_reads_before_filtering

if not merged_wide.total_nodups.equals(merged_wide.cis + merged_wide.trans):
    raise ValueError("total_nodups")
if not merged_wide.total_mapped.equals(merged_wide.total_dups + merged_wide.total_nodups):
    raise ValueError("total_mapped")
if not merged_wide.total.equals(merged_wide.total_mapped + merged_wide.total_single_sided_mapped + merged_wide.total_unmapped):
    raise ValueError("total")

merged_wide.to_csv("pairs_stats.wide.csv")
merged_wide[["total_before_filtering","total","total_mapped","total_nodups","cis"]].to_csv("pairs_stats.wide.subset1.csv")

plt.figure()
merged_wide[["cis", "trans"]].plot.bar(stacked=True)
plt.subplots_adjust(bottom=0.5)
plt.savefig("pairs_stats.cis_trans.svg")

plt.figure()
merged_wide["cis_20kb+"].plot.bar()
plt.subplots_adjust(bottom=0.5)
plt.savefig("pairs_stats.cis20kb.svg")
