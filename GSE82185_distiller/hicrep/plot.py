#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
samples = [
    "GSE82185_early2C_No1",
    "GSE82185_early2C_No2",
    "GSE82185_early2C_No3",
    "GSE82185_late2C_No1",
    "GSE82185_late2C_No2",
    "GSE82185_late2C_No3",
    "GSE82185_late2C_No4",
    ]
sample_names = [
    "GSE82185 Early 2-cell #1",
    "GSE82185 Early 2-cell #2",
    "GSE82185 Early 2-cell #3",
    "GSE82185 Late 2-cell #1",
    "GSE82185 Late 2-cell #2",
    "GSE82185 Late 2-cell #3",
    "GSE82185 Late 2-cell #4",
    ]
similarities = np.empty((len(samples), len(samples)))
for i in range(len(samples)):
    for j in range(i, len(samples)):
        hicrep_path = f"hicrep_out.{samples[i]}.{samples[j]}.csv"
        d = pd.read_csv(hicrep_path)
        similarities[i, j] = d.iloc[:,1].median()
        similarities[j, i] = similarities[i, j]
print(similarities)
fig, ax = plt.subplots()
#image = ax.matshow(similarities, cmap="inferno")
image = ax.matshow(similarities, cmap="coolwarm", vmin=-1, vmax=1)
ax.set_xticks(ticks=range(len(samples)), labels=sample_names, rotation=90)
ax.set_yticks(ticks=range(len(samples)), labels=sample_names, rotation=0)
fig.colorbar(image, label="Stratum-adjusted Correlation Coefficient")
for ((x, y), v) in np.ndenumerate(similarities):
    ax.text(y, x, f"{v:.2f}", ha="center", va="center")
fig.tight_layout()
fig.savefig("plot.matrix.png", dpi=300)
