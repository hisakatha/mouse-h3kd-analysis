#!/usr/bin/env python3
"""Get loop strengths
This script gets loop strengths from .mcool files

References
----------
Flyamer2017: Flyamer et al. (2017). Single-nucleus Hi-C reveals unique chromatin
reorganization at oocyte-to-zygote transition. Nature, doi:10.1038/nature21711

https://coolpuppy.readthedocs.io/en/v1.0.0/Examples/TAD_score.html
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from coolpuppy import coolpup
from coolpuppy.lib.numutils import get_domain_score
from coolpuppy.lib.puputils import accumulate_values
from coolpuppy import plotpup
import cooler
import bioframe
import cooltools
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

def get_loop_strength(mcool_path, resolution, chrom_arms, loops):
    clr = cooler.Cooler(f"{mcool_path}::resolutions/{resolution}")
    #expected = expected_cis(clr, ignore_diags=0, view_df=chrom_arms, chunksize=1000000)
    cc = coolpup.CoordCreator(loops, resolution=resolution, features_format='bedpe')
    pu = coolpup.PileUpper(clr, cc, view_df=chrom_arms, control=True)
    #pup = pu.pileupsWithControl(postprocess_func=add_domain_score,
    #                        extra_sum_funcs={'domain_score': extra_sum_func})
    pup = pu.pileupsWithControl()
    return pup

"""
genome_id = "mm10"
chromsizes = bioframe.fetch_chromsizes(genome_id)
chrom_cens = bioframe.fetch_centromeres(genome_id)
chrom_arms = bioframe.make_chromarms(chromsizes, chrom_cens)
chrom_arms = bioframe.make_viewframe(chrom_arms)
"""
chrom_arms = None

#bedpe_header = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
bedpe_header = dict(enumerate(["chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2"]))
loops = pd.read_table("./GSE63525_mouse_lymphoblasts_HiCCUPS_looplist.txt.mm10.bedpe", header=None)
loops.rename(columns=bedpe_header, inplace=True)
if not all(loops.chrom1 == loops.chrom2):
    raise ValueError("Input loop positions should satisfy chrom1 == chrom2")

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

loops["chrom1"] = loops["chrom1"].map(ensembl_mm10_from_ucsc)
loops["chrom2"] = loops["chrom2"].map(ensembl_mm10_from_ucsc)

resolution = 10000
#pileup_size = 99

clr_path_2c_17hpi = "./coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool"
clr_path_2c_control = "./coolers_library_group/2c_control.mm10.no_filter.1000.mcool"
clr_path_2c_h3kd = "./coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool"
clr_paths = {
    "2C 17hpi": clr_path_2c_17hpi,
    "2C Control": clr_path_2c_control,
    "2C H3.1/3.2 KD": clr_path_2c_h3kd,
}

#clr_path_oocyte = "./oocyte_SN_3.mm10.no_filter.1000.mcool"
#clr_paths["Oocyte SN #3"] = clr_path_oocyte

#clr_path_gse82185_2c = "./GSE82185_late2C.mm10.no_filter.1000.mcool"
#clr_paths["2C GSE82185"] = clr_path_gse82185_2c

#clr_path_gse100569_1c_g2 = "./GSE100569_G2.mm10.no_filter.1000.mcool"
#clr_paths["1C G2 GSE100569"] = clr_path_gse100569_1c_g2

samples = list(clr_paths.keys())

height = 3

pup_dict = {k: get_loop_strength(v, resolution, chrom_arms, loops) for (k,v) in clr_paths.items()}
pups = pd.concat([v.assign(Sample=k) for (k,v) in pup_dict.items()])

plotpup.plot(pups, height=height, cols="Sample", col_order=samples, score="", plot_ticks=True)
plt.savefig("loop_strength.heatmap.svg")

exit(0)
############ END ############

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
print(test_score_diff(pups, "2C Control", "2C H3.1/3.2 KD"))
print("Rank sum test: Control vs 17hpi")
print(test_score_diff(pups, "2C Control", "2C 17hpi"))
print("Rank sum test: 17hpi vs KD")
print(test_score_diff(pups, "2C 17hpi", "2C H3.1/3.2 KD"))

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
plt.savefig("tad_strength.ecdf.svg")
