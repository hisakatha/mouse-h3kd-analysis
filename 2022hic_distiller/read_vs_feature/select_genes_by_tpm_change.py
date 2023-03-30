#!/usr/bin/env python3
import pandas as pd
import numpy as np

tpms = pd.read_table("Embryonicgene_5stage_average_tpm.tsv").set_index("Geneid")
#clusters = pd.concat([
#    pd.read_table("tpm_3stage_zygoticgene_Noinjection_up.tsv", usecols=[0]).assign(cluster="Up"),
#    pd.read_table("tpm_3stage_zygoticgene_Noinjection_nochange.tsv", usecols=[0]).assign(cluster="NoChange"),
#    pd.read_table("tpm_3stage_zygoticgene_Noinjection_down.tsv", usecols=[0]).assign(cluster="Down"),
#]).set_index("Geneid")
#clusters.join(tpms, how="inner", validate="one_to_one")

clusters2 = pd.concat([
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_up.tsv").assign(cluster="Up"),
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_nochange.tsv").assign(cluster="NoChange"),
    pd.read_table("tpm_3stage_zygoticgene_Noinjection_down.tsv").assign(cluster="Down"),
]).set_index("Geneid")
genes = clusters2[["cluster", "Mll"]].join(tpms, how="inner", validate="one_to_one")
genes["log2fc_2c_18hpi_vs_Mll"] = np.log2(genes["2c_18hpi"] / genes["Mll"])
genes["log2fc_2c_Control_vs_2c_18hpi"] = np.log2(genes["2c_Control"] / genes["2c_18hpi"])
genes["log2fc_2c_H3_1_3_2KD_vs_2c_Control"] = np.log2(genes["2c_H3_1_3_2KD"] / genes["2c_Control"])

cluster_up_target_flags = (genes["cluster"] == "Up") & np.isfinite(genes["log2fc_2c_18hpi_vs_Mll"]) & \
    (np.abs(genes["log2fc_2c_18hpi_vs_Mll"]) < 1) & (genes["log2fc_2c_H3_1_3_2KD_vs_2c_Control"] < 0)
cluster_up_targets = genes[cluster_up_target_flags].sort_values(by="log2fc_2c_Control_vs_2c_18hpi", ascending=False)[:10].index.to_list()
cluster_down_target_flags = (genes["cluster"] == "Down") & np.isfinite(genes["log2fc_2c_H3_1_3_2KD_vs_2c_Control"])
cluster_down_targets = genes[cluster_down_target_flags].sort_values(by="log2fc_2c_H3_1_3_2KD_vs_2c_Control", ascending=False)[:10].index.to_list()

bed_cols = ["chr", "start", "end", "name", "score", "strand"]
gene_regions = pd.read_table("promoter.embryonic_igenomes_ucsc.ensembl.bed", names=bed_cols)
#gene_regions[gene_regions["name"].duplicated()]

genes2 = genes.join(gene_regions.set_index("name"), how="outer")
genes2.loc[cluster_up_targets].to_csv("select_genes_by_tpm_change.cluster_up.csv", index_label="name")
genes2.loc[cluster_down_targets].to_csv("select_genes_by_tpm_change.cluster_down.csv", index_label="name")
