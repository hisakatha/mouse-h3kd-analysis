#!/usr/bin/env python3
import pandas as pd
import pickle
import gzip
import matplotlib.pyplot as plt

tpms = pd.read_table("Embryonicgene_5stage_average_tpm.tsv").set_index("Geneid")

prefix = "get_enhancer_ratio_intrachr_100k"
with gzip.open(prefix + ".pkl.gz", "rb") as fp:
    feature1_summaries = pickle.load(fp)

hic_name_to_tpm_name = {
    "2C 17hpi 1001bp embryonic genes": "2c_18hpi",
    "2C control 1001bp embryonic genes": "2c_Control",
    "2C H3.1/3.2KD 1001bp embryonic genes": "2c_H3_1_3_2KD",
    "2C 17hpi 2001bp embryonic genes": "2c_18hpi",
    "2C control 2001bp embryonic genes": "2c_Control",
    "2C H3.1/3.2KD 2001bp embryonic genes": "2c_H3_1_3_2KD",
}

fig, axs = plt.subplots(len(feature1_summaries), 1, squeeze=False, sharex=True, sharey=True, figsize=(3.5,11))
for (idx, (sample, feature1_summary)) in enumerate(feature1_summaries.items()):
    ax = axs[idx, 0]
    d2 = feature1_summary.join(tpms)
    sample_name_in_tpm = hic_name_to_tpm_name[sample]
    #ax.scatter(d2["2c_Control"], d2["num_reads_with_feature2"], alpha=0.1)
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    d2 = d2[(d2[sample_name_in_tpm] > 0) & (d2["num_reads_with_feature2"] > 0)]
    hb = ax.hexbin(d2[sample_name_in_tpm], d2["num_reads_with_feature2"], xscale="log", yscale="linear", gridsize=25, extent=(0,4,10,60), bins="log")
    #hb = ax.hexbin(d2[sample_name_in_tpm], d2["num_reads_with_feature2"], xscale="log", yscale="log", gridsize=25, bins="log")
    ax.set_xticks([1,10**2,10**4])
    ax.set_title(sample, fontsize="medium")
    ax.set_box_aspect(1)
    fig.colorbar(hb, ax=ax)

fig.supxlabel("TPM")
fig.supylabel("#promoter-enhancer reads")
fig.savefig(prefix + ".plot_interaction_vs_expr.png", dpi=300)
