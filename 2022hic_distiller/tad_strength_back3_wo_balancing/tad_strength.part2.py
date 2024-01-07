#!/usr/bin/env python3
"""Plot figures from the outputs of tad_strength.py
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import seaborn as sns

# TAD strengths calculated with tad_strength.py
pups_path = "tad_strength.pandas_pkl.gz"
pups = pd.read_pickle(pups_path)
tad_strength_early = pups.loc[pups.Sample == "Early 2-cell"].domain_score[0]
tad_strength_control = pups.loc[pups.Sample == "Late 2-cell Control"].domain_score[0]
tad_strength_h3kd = pups.loc[pups.Sample == "Late 2-cell H3.1/3.2 KD"].domain_score[0]

if len(tad_strength_control) != len(tad_strength_h3kd):
    raise ValueError("The number of TAD strength must be same")
if len(tad_strength_control) != len(tad_strength_early):
    raise ValueError("The number of TAD strength must be same")
if not all(np.isfinite(tad_strength_control)):
    raise ValueError("All of the TAD strengths must be infinite")
if not all(np.isfinite(tad_strength_h3kd)):
    raise ValueError("All of the TAD strengths must be infinite")

##### Get correlations between samples #####
tad_pearsonr = stats.pearsonr(np.log10(tad_strength_control), np.log10(tad_strength_h3kd))
print(tad_pearsonr)
with open("tad_strength.cor_h3kd.txt", "w") as f:
    f.write(f"Pearson's correlation test: Control vs. H3KD\nPearson's correlation coefficient: {tad_pearsonr.statistic}\nP-value (two-sided): {tad_pearsonr.pvalue}\n")

tad_strength_table = pd.DataFrame({"log10(Late 2-cell Control)": np.log10(tad_strength_control),
                                   "log10(Late 2-cell H3KD)": np.log10(tad_strength_h3kd)})

fig, axs = plt.subplots(squeeze=False)
ax = axs[0, 0]
sns.regplot(x="log10(Late 2-cell Control)", y="log10(Late 2-cell H3KD)", data=tad_strength_table,
            scatter_kws={"alpha":0.2}, ax=ax)
ax.axline((0, 0), slope=1, color="black", linestyle="--")
ax.set_xlabel(r"$\log_{10}$(TAD strengths of Late 2-cell Control)")
ax.set_ylabel(r"$\log_{10}$(TAD strengths of Late 2-cell H3KD)")
fig.savefig("tad_strength.cor_h3kd_scatter.svg")

fig, axs = plt.subplots(squeeze=False)
ax = axs[0, 0]
ax.set_aspect("equal", "box")
#ax.axis("equal")
hb = ax.hexbin(x="log10(Late 2-cell Control)", y="log10(Late 2-cell H3KD)", data=tad_strength_table,
          gridsize=33, edgecolors="none", bins=None, mincnt=1)
cbar = fig.colorbar(hb, label="Counts", shrink=0.8)
#cbar.set_ticks(np.unique(np.concatenate((cbar.get_ticks(), [1]))).tolist())
cbar.set_ticks([1,10,20,30,40,50])
sns.regplot(x="log10(Late 2-cell Control)", y="log10(Late 2-cell H3KD)", data=tad_strength_table,
            scatter=False, ax=ax, color="red")
ax.axline((0, 0), slope=1, color="black", linestyle="--")
ax.set_xlabel(r"$\log_{10}$(TAD strengths of Late 2-cell Control)")
ax.set_ylabel(r"$\log_{10}$(TAD strengths of Late 2-cell H3KD)")
fig.savefig("tad_strength.cor_h3kd_hex.svg")

##### Get differences between samples #####
log10_tad_strength_ratios = np.log10(np.array(tad_strength_h3kd) / np.array(tad_strength_control))
tad_signed_rank = stats.wilcoxon(log10_tad_strength_ratios)
print(tad_signed_rank)
with open("tad_strength.ratio_h3kd_signed_rank.txt", "w") as f:
    f.write(f"Wilcoxon signed rank test: Control vs. H3KD\nStatistic (two-sided): {tad_signed_rank.statistic}\nP-value (two-sided): {tad_signed_rank.pvalue}\n")

fig, axs = plt.subplots(squeeze=False)
ax = axs[0, 0]
sns.boxplot(log10_tad_strength_ratios)
ax.axhline(0, color="black", linestyle="--")
ax.set_box_aspect(4)
ax.set_xticks([])
ax.set_ylabel(r"$\log_{10}$(TAD str. of H3KD) - $\log_{10}$(TAD str. of Control)")
fig.savefig("tad_strength.ratio_h3kd_box.svg")

fig, axs = plt.subplots(squeeze=False)
ax = axs[0, 0]
sns.histplot(log10_tad_strength_ratios, bins=np.arange(-1.6,1.6,0.1))
ax.set_xlabel(r"$\log_{10}$(TAD str. of H3KD) - $\log_{10}$(TAD str. of Control)")
fig.savefig("tad_strength.ratio_h3kd_hist.svg")
