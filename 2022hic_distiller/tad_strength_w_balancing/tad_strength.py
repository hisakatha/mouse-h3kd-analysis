#!/usr/bin/env python3
"""Get TAD strengths
This script gets TAD strengths from .mcool files

References
----------
Flyamer2017: Flyamer et al. (2017). Single-nucleus Hi-C reveals unique chromatin
reorganization at oocyte-to-zygote transition. Nature, doi:10.1038/nature21711

https://coolpuppy.readthedocs.io/en/v1.0.0/Examples/TAD_score.html
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from coolpuppy import coolpup
from coolpuppy.lib.numutils import get_domain_score
from coolpuppy.lib.puputils import accumulate_values
from coolpuppy import plotpup
import cooler
#import bioframe
#import cooltools
#from cooltools.lib import io
from cooltools import expected_cis
#from cooltools.lib import plotting
from itertools import compress
from scipy import stats
import seaborn as sns

def add_domain_score(snippet):
    # Calculates domain score for each snippet according to Flyamer et al., 2017
    snippet['domain_score'] = get_domain_score(snippet['data'])
    return snippet

def extra_sum_func(dict1, dict2):
    return accumulate_values(dict1, dict2, 'domain_score')

def get_tad_strength(mcool_path, resolution, pileup_size, chrom_arms, tads):
    clr = cooler.Cooler(f"{mcool_path}::resolutions/{resolution}")
    expected = expected_cis(clr, ignore_diags=0, view_df=chrom_arms, chunksize=1000000)
    cc = coolpup.CoordCreator(tads, resolution=resolution, features_format='bed', local=True, rescale_flank=1)
    pu = coolpup.PileUpper(clr, cc, expected=expected, view_df=chrom_arms, ignore_diags=0, rescale_size=pileup_size, rescale=True)
    pup = pu.pileupsWithControl(postprocess_func=add_domain_score,
                            extra_sum_funcs={'domain_score': extra_sum_func})
    return pup

#genome_id = "mm10"
#chromsizes = bioframe.fetch_chromsizes(genome_id)
#chrom_cens = bioframe.fetch_centromeres(genome_id)
#chrom_arms = bioframe.make_chromarms(chromsizes, chrom_cens)
#chrom_arms = bioframe.make_viewframe(chrom_arms)
chrom_arms = None

tads = pd.read_table("../GSE63525_mouse_lymphoblasts_Arrowhead_domainlist.txt.mm10.bed", names=["chrom", "start", "end"])

ensembl_mm10_from_ucsc = {
    "chr1": "1",
    "chr2": "2",
    "chr3": "3",
    "chr4": "4",
    "chr5": "5",
    "chr6": "6",
    "chr7": "7",
    "chr8": "8",
    "chr9": "9",
    "chr10": "10",
    "chr11": "11",
    "chr12": "12",
    "chr13": "13",
    "chr14": "14",
    "chr15": "15",
    "chr16": "16",
    "chr17": "17",
    "chr18": "18",
    "chr19": "19",
    "chrM": "MT",
    "chrX": "X",
    "chrY": "Y",
}

tads["chrom"] = tads["chrom"].map(ensembl_mm10_from_ucsc)

resolution = 10000
pileup_size = 99

clr_path_2c_17hpi = "../coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool"
clr_path_2c_control = "../coolers_library_group/2c_control.mm10.no_filter.1000.mcool"
clr_path_2c_h3kd = "../coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool"
clr_paths = {
    "Early 2-cell": clr_path_2c_17hpi,
    "Late 2-cell Control": clr_path_2c_control,
    "Late 2-cell H3.1/3.2 KD": clr_path_2c_h3kd,
}

#clr_path_oocyte = "./oocyte_SN_3.mm10.no_filter.1000.mcool"
#clr_paths["Oocyte SN #3"] = clr_path_oocyte

#clr_path_gse82185_2c = "./GSE82185_late2C.mm10.no_filter.1000.mcool"
#clr_paths["2C GSE82185"] = clr_path_gse82185_2c
samples = list(clr_paths.keys())

height = 3

pups_path = "tad_strength.pandas_pkl.gz"
if not os.path.exists(pups_path):
    pup_dict = {k: get_tad_strength(v, resolution, pileup_size, chrom_arms, tads) for (k,v) in clr_paths.items()}
    pups = pd.concat([v.assign(Sample=k) for (k,v) in pup_dict.items()])
    pups.to_pickle(pups_path)
else:
    pups = pd.read_pickle(pups_path)

plotpup.plot(pups, height=height, cols="Sample", col_order=samples, score="", plot_ticks=False)
plt.annotate("Enrichment", xy=(0.999, 0.5), xycoords="figure fraction", rotation=270, ha="right", va="center")
plt.savefig("tad_strength.heatmap.svg")

plotpup.plot(pups, height=height, cols="Sample", col_order=samples, score="", plot_ticks=False, vmax=2)
plt.annotate("Enrichment", xy=(0.999, 0.5), xycoords="figure fraction", rotation=270, ha="right", va="center")
plt.savefig("tad_strength.heatmap_vmax2.svg")

plotpup.plot(pups, height=height, cols="Sample", col_order=samples, score="", plot_ticks=False, vmax=1.5)
plt.annotate("Enrichment", xy=(0.999, 0.5), xycoords="figure fraction", rotation=270, ha="right", va="center")
plt.savefig("tad_strength.heatmap_vmax1.5.svg")

pups["domain_score_finite"] = pups["domain_score"].apply(lambda x: list(compress(x, np.isfinite(x))))
pups["domain_score_finite_mean"] = pups["domain_score_finite"].apply(np.mean)

fig, axs = plt.subplots(len(pups), 1, sharex=True, sharey=False)
for i in range(len(pups)):
    sample = samples[i]
    print(sample)
    if sum(pups["Sample"] == sample) > 1:
        raise ValueError("Unexpected sample number")
    sample_pup = pups.loc[pups["Sample"] == sample, :]
    print("max domain_score", max(sample_pup["domain_score"][0]))
    print("max domain_score_finite", max(sample_pup["domain_score_finite"][0]))
    axs[i].hist(sample_pup["domain_score_finite"], bins=np.linspace(0,4,51))
    axs[i].axvline(sample_pup["domain_score_finite_mean"][0], color="red")
    axs[i].set_title(sample, loc="right", y=0.6)
    if i in range(3):
        axs[i].set_ylim(0, 260)
    axs[i].set_xlim(0, 4)
fig.supylabel("Count")
fig.supxlabel("TAD strength")
fig.suptitle("Distribution of TAD strengths")
fig.savefig("tad_strength.hist.svg")

plt.figure()
plt.boxplot(pups["domain_score_finite"], labels=pups["Sample"])
plt.ylim(-0.5, 10)
plt.ylabel("TAD strength")
plt.savefig("tad_strength.boxplot.svg")
#dict1 = {row.Sample:list(compress(row.domain_score, np.isfinite(row.domain_score))) for row in pups2.itertuples()}

def test_score_diff(data, sample1, sample2, score="domain_score_finite"):
    return stats.mannwhitneyu(data.loc[data["Sample"] == sample1, score][0], data.loc[data["Sample"] == sample2, score][0])
print("Rank sum test: Control vs KD")
print(test_score_diff(pups, "Late 2-cell Control", "Late 2-cell H3.1/3.2 KD"))
print("Rank sum test: Control vs 17hpi")
print(test_score_diff(pups, "Late 2-cell Control", "Early 2-cell"))
print("Rank sum test: 17hpi vs KD")
print(test_score_diff(pups, "Early 2-cell", "Late 2-cell H3.1/3.2 KD"))

plt.figure()
#pups.plot.bar(x="Sample", y="domain_score_finite_mean", legend=False)
plt.bar(x=pups["Sample"], height=pups["domain_score_finite_mean"])
plt.xticks(rotation=15)
plt.xlabel("Sample")
plt.ylabel("Average TAD strength")
plt.savefig("tad_strength.bar.svg", bbox_inches='tight', pad_inches=0.1)

pups_long = pd.concat([pd.DataFrame({
    "domain_score_finite": pups.iloc[i]["domain_score_finite"],
    "Sample": pups.iloc[i]["Sample"]})
    for i in range(len(pups))])

plt.figure()
#sns.swarmplot(pups_long, x="Sample", y="domain_score_finite") # This is very slow!
sns.stripplot(pups_long, x="Sample", y="domain_score_finite", alpha=0.2, jitter=0.3)
plt.ylim(0, 4)
plt.ylabel("TAD strength")
plt.savefig("tad_strength.jitter.png")

plt.figure()
sns.ecdfplot(pups_long, hue="Sample", x="domain_score_finite", linestyle="-")
plt.xlabel("TAD strength")
plt.xlim(0, 4)
#plt.ylabel("Cumulative probability")
plt.savefig("tad_strength.ecdf.svg")

# Proportion of TAD strength larger than a threshold
from collections import namedtuple
domain_score_thres = 2.5
TAD_Count = namedtuple("TAD_Count", ["Sample", "total", "num_ge_thres"])
tad_counts = pd.DataFrame([TAD_Count(pups.iloc[i]["Sample"], len(pups.iloc[i]["domain_score_finite"]),
                                     sum(score >= domain_score_thres for score in pups.iloc[i]["domain_score_finite"])) for i in range(len(pups))])
tad_counts["num_lt_thres"] = tad_counts["total"] - tad_counts["num_ge_thres"]

def compare_tad_counts(tad_counts, sample1, sample2):
    tad_counts_sub = tad_counts.loc[tad_counts["Sample"].isin([sample1, sample2]),["num_lt_thres","num_ge_thres"]]
    #stats.contingency.chi2_contingency(tad_counts_sub)
    (_, p) = stats.fisher_exact(tad_counts_sub)
    print(f"{sample1},{sample2},{p}")

print(f"Category size diff for TAD strength threshold of {domain_score_thres}")
compare_tad_counts(tad_counts, "Late 2-cell Control", "Early 2-cell")
compare_tad_counts(tad_counts, "Late 2-cell Control", "Late 2-cell H3.1/3.2 KD")

plt.figure()
tad_counts["frac_lt_thres"] = tad_counts["num_lt_thres"] / tad_counts["total"]
tad_counts["frac_ge_thres"] = tad_counts["num_ge_thres"] / tad_counts["total"]
#tad_counts.plot.bar(x="Sample", y=["frac_lt_thres", "frac_ge_thres"], ylabel="Proportion", stacked=True, rot=0)
#plt.legend([f"TAD strength < {domain_score_thres}", f"TAD strength ≧ {domain_score_thres}"], loc="lower left")
tad_counts.plot.bar(x="Sample", y=["frac_ge_thres"], ylabel=f"Proportion of TAD strengths ≧ {domain_score_thres}", stacked=True, rot=0, legend=False)
plt.savefig("tad_strength.proportion_w_thres.svg")

plt.figure()
tad_counts.plot.bar(x="Sample", y=["num_ge_thres"], ylabel=f"The number of TAD strengths ≧ {domain_score_thres}", stacked=True, rot=0, legend=False)
plt.savefig("tad_strength.num_ge_thres.svg")
