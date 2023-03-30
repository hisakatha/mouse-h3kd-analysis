#!/usr/bin/env python3
import os
import gzip
import math
import pickle
import pandas as pd
# coolbox requires 'bgzip' and 'tabix' to be available in PATH
import coolbox
import coolbox.api as cb
#import matplotlib.pyplot as plt
import svgutils.compose as svg

sample_names = {
        "2c_17hpi_No1": "Early 2-cell #1",
        "2c_17hpi_No2": "Early 2-cell #2",
        "2c_control_No1": "Late 2-cell Control #1",
        "2c_control_No2": "Late 2-cell Control #2",
        "2c_H3-1_and_3-2KD_No1": "Late 2-cell H3.1/3.2 KD #1",
        "2c_H3-1_and_3-2KD_No2": "Late 2-cell H3.1/3.2 KD #2",
        "2c_17hpi": "Early 2-cell",
        "2c_control": "Late 2-cell Control",
        "2c_H3-1_and_3-2KD": "Late 2-cell H3.1/3.2 KD",
        }

sample_paths = {
        "2c_17hpi_No1": "../results/coolers_library/2c_17hpi_No1.mm10.no_filter.1000.mcool",
        "2c_17hpi_No2": "../results/coolers_library/2c_17hpi_No2.mm10.no_filter.1000.mcool",
        "2c_control_No1": "../results/coolers_library/2c_control_No1.mm10.no_filter.1000.mcool",
        "2c_control_No2": "../results/coolers_library/2c_control_No2.mm10.no_filter.1000.mcool",
        "2c_H3-1_and_3-2KD_No1": "../results/coolers_library/2c_H3-1_and_3-2KD_No1.mm10.no_filter.1000.mcool",
        "2c_H3-1_and_3-2KD_No2": "../results/coolers_library/2c_H3-1_and_3-2KD_No2.mm10.no_filter.1000.mcool",
        "2c_17hpi": "../results/coolers_library_group/2c_17hpi.mm10.no_filter.1000.mcool",
        "2c_control": "../results/coolers_library_group/2c_control.mm10.no_filter.1000.mcool",
        "2c_H3-1_and_3-2KD": "../results/coolers_library_group/2c_H3-1_and_3-2KD.mm10.no_filter.1000.mcool",
        }

#sample_ids = list(sample_names.keys())
replicate_sample_ids = [
        "2c_17hpi_No1",
        "2c_17hpi_No2",
        "2c_control_No1",
        "2c_control_No2",
        "2c_H3-1_and_3-2KD_No1",
        "2c_H3-1_and_3-2KD_No2",
        ]
merged_sample_ids = [
        "2c_17hpi",
        "2c_control",
        "2c_H3-1_and_3-2KD",
        ]

genes_dtypes = {
        "chr": str,
        }
genes_cluster_up = pd.read_csv("select_genes_by_tpm_change.cluster_up.csv", dtype=genes_dtypes)
genes_cluster_down = pd.read_csv("select_genes_by_tpm_change.cluster_down.csv", dtype=genes_dtypes)
bed_display = "stacked"
cool_vmin1 = 0

# For resolution = 100k
cool_vmax1 = 12
resolution = 100000
total_plot_length = 10_000_000

# For resolution = 25k
#cool_vmax1 = 11
#resolution = 25000
#total_plot_length = 2_000_000

# annotations compressed with bgzip
promoter_path = "promoter.embryonic_igenomes_ucsc.ensembl.bed.bgz"
bed6_cols = ["chr", "start", "end", "name", "score", "strand"]
promoters = pd.read_table(promoter_path, compression="gzip", header=None, names=bed6_cols)
enhancer_path = "ATAC_early2C_only_DMSO.extend.ensembl.bed.bgz"
chrom_lengths_table = pd.read_table("Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai", usecols=[0,1], names=["chr", "length"]).set_index("chr")
chrom_lengths = chrom_lengths_table.to_dict()["length"]

with gzip.open("get_enhancer_ratio_intrachr_1k_replicate.pkl.gz", "rb") as fp:
    promoter_summaries = pickle.load(fp)
with gzip.open("get_enhancer_ratio_intrachr_1k.pkl.gz", "rb") as fp:
    promoter_summaries2 = pickle.load(fp)
promoter_summaries.update(promoter_summaries2)

promoter_summaries_keys = {
        "2c_17hpi_No1": "2C 17hpi #1 2001bp embryonic genes",
        "2c_17hpi_No2": "2C 17hpi #2 2001bp embryonic genes",
        "2c_control_No1": "2C control #1 2001bp embryonic genes",
        "2c_control_No2": "2C control #2 2001bp embryonic genes",
        "2c_H3-1_and_3-2KD_No1": "2C H3.1/3.2KD #1 2001bp embryonic genes",
        "2c_H3-1_and_3-2KD_No2": "2C H3.1/3.2KD #2 2001bp embryonic genes",
        "2c_17hpi": "2C 17hpi 2001bp embryonic genes",
        "2c_control": "2C control 2001bp embryonic genes",
        "2c_H3-1_and_3-2KD": "2C H3.1/3.2KD 2001bp embryonic genes",
        }

def plot_triangular(cool_path, promoter_path, enhancer_path, region1, sample_name_fig, sample_name_out, region_name_fig, region_name_out,
    balance=False, bed_display="stacked", promoter_highlight=None, cool_vmin="auto", cool_vmax="auto"):
    bed_color = "blue"
    frame_promoter = cb.Spacer(1) + cb.BED(promoter_path, color=bed_color, display=bed_display, labels=False) + cb.Title("Promoters")
    if promoter_highlight is not None:
        #frame_promoter += cb.HighLights(promoter_highlight)
        frame_promoter += cb.Vlines(promoter_highlight)
    frame_enhancer = cb.Spacer(1) + cb.BED(enhancer_path, color=bed_color, display=bed_display, labels=False) + cb.Title("Enhancers")
    frame = cb.Spacer(0.5) + cb.XAxis() + cb.Cool(cool_path, resolution=resolution, style="triangular",
        balance=balance, transform="log2", cmap="plasma", min_value=cool_vmin, max_value=cool_vmax)
    frame += frame_promoter + frame_enhancer + cb.FrameTitle(f"{region_name_fig} ({region1}): {sample_name_fig}")
    fig = frame.plot(region1)
    #print(f"plot_triangular: {type(fig)}")
    fig.savefig(f"coolbox_plot.triangular.{region_name_out}.{sample_name_out}.png")
    return fig

def plot_2d(cool_path, promoter_path, enhancer_path, region1, region2, sample_name_fig, sample_name_out, region_name_fig, region_name_out,
    balance=False, bed_display="stacked", promoter_highlight=None, cool_vmin="auto", cool_vmax="auto"):
    bed_color = "blue"
    if (region1 is None) and (region2 is None):
        raise ValueError("Either region1 or region2 must have a region string")
    elif region1 is None:
        region1 = region2
        region1_label = "Promoters"
        region2_label = f"Enhancers ({region2})"
    elif region2 is None:
        region2 = region1
        region1_label = f"Promoters ({region1})"
        region2_label = "Enhancers"
    else:
        region1_label = f"Promoters ({region1})"
        region2_label = f"Enhancers ({region2})"
    frame_promoter = cb.Spacer(0.5) + cb.XAxis() + cb.BED(promoter_path, color=bed_color, display=bed_display, labels=False) +\
        cb.FrameTitle(f"{region_name_fig}: {sample_name_fig}\n{region1_label}")
    if promoter_highlight is not None:
        #frame_promoter += cb.HighLights(promoter_highlight)
        frame_promoter += cb.Vlines(promoter_highlight)
    frame_enhancer = cb.XAxis() + cb.BED(enhancer_path, color=bed_color, display=bed_display, labels=False) + cb.FrameTitle(region2_label)
    sub_frames = {
        "top": frame_promoter,
        "left": frame_enhancer,
        "bottom": cb.Spacer(2),
        "right": cb.Spacer(2),
        }
    frame_center = cb.Cool(cool_path, resolution=resolution, style="matrix", color_bar="no",
        min_value=cool_vmin, max_value=cool_vmax, balance=balance, transform="log2", cmap="plasma")
    #frame_center = cb.Cool(cool_path, resolution=resolution, style="matrix", color_bar="no", balance=balance, transform="log2", cmap="plasma")
    frame = cb.JointView(frame_center, **sub_frames, space=0, padding_left=0.5)
    fig = frame.plot(region1, region2)
    #print(f"plot_2d: {type(fig)}")
    #print(f"plot_2d: {fig.height}")
    # This figure object is from svgutils package and can be saved only in svg format
    fig.save(f"coolbox_plot.2d.{region_name_out}.{sample_name_out}.svg")
    return fig

def get_extended_region(gene_info, total_length, chrom_lengths):
    mid = (gene_info.start + gene_info.end) // 2
    extended_start = mid - (total_length // 2)
    extended_start = extended_start if (extended_start >= 0) else 0
    extended_end = mid + ((total_length + 1) // 2) - 1
    extended_end = extended_end if (extended_end <= chrom_lengths[gene_info.chr]) else chrom_lengths[gene_info.chr]
    return f"{gene_info.chr}:{extended_start}-{extended_end}"


def plot_all_samples(region, region_name_fig, region_name_out, sample_ids, all_name_out):
    region1_gene = [(region.chr, (region.start + region.end)//2)]
    region1_plot = get_extended_region(region, total_plot_length, chrom_lengths)
    focused_promoter_path = f"coolbox_plot.promoter.{region_name_out}.bed"
    focused_compressed_promoter_path = f"{focused_promoter_path}.bgz"
    promoters[promoters["name"] == region.name].to_csv(focused_promoter_path, sep="\t", header=False, index=False)
    os.system(f"bgzip -c {focused_promoter_path} > {focused_compressed_promoter_path} && tabix -p bed {focused_compressed_promoter_path}")
    figures_triangular = []
    figures_2d = []
    for sample_id in sample_ids:
        # enhancer-promoter (EP) reads
        ep_read_info = promoter_summaries[promoter_summaries_keys[sample_id]].loc[region.name]
        ep_read_num_raw = ep_read_info["num_reads_with_feature2"]
        ep_read_num = int(ep_read_num_raw)
        if not math.isclose(ep_read_num_raw, ep_read_num):
            raise ValueError("num_reads_with_feature2 must be (close to) an integer")
        ep_read_ratio = ep_read_info["ratio"]
        region_name_fig_for_sample = f"{region_name_fig} (#EP = {ep_read_num:d}, #EP/#P = {ep_read_ratio:.3g})"
        sample_path = sample_paths[sample_id]
        sample_name = sample_names[sample_id]
        fig_triangular = plot_triangular(sample_path, focused_compressed_promoter_path, enhancer_path,
            region1_plot, sample_name, sample_id, region_name_fig_for_sample, region_name_out,
            bed_display=bed_display, promoter_highlight=None)
        figures_triangular.append(svg.MplFigure(fig_triangular).scale(0.8))
        fig_2d = plot_2d(sample_path, focused_compressed_promoter_path, enhancer_path,
            region1_plot, None, sample_name, sample_id, region_name_fig_for_sample, region_name_out,
            bed_display=bed_display, promoter_highlight=None, cool_vmin=cool_vmin1, cool_vmax=cool_vmax1)
        figures_2d.append(fig_2d)
    # figure width in svgutils Unit object
    single_image_width = figures_2d[0].width
    svg.Figure(str(single_image_width * len(sample_ids)), str(single_image_width), *figures_triangular).tile(len(sample_ids), 1)\
        .save(f"coolbox_plot.triangular.{region_name_out}.{all_name_out}.svg")
    svg.Figure(str(single_image_width * len(sample_ids)), str(single_image_width), *figures_2d).tile(len(sample_ids), 1)\
        .save(f"coolbox_plot.2d.{region_name_out}.{all_name_out}.svg")

for (i, gene) in enumerate(genes_cluster_up.itertuples(), start=1):
    region_name_fig = f"Cluster UP rank {i}: {gene.name}"
    region_name_out = f"up_rank{i:03}"
    plot_all_samples(gene, region_name_fig, region_name_out, replicate_sample_ids, "all_replicate_samples")
    plot_all_samples(gene, region_name_fig, region_name_out, merged_sample_ids, "all_merged_samples")
for (i, gene) in enumerate(genes_cluster_down.itertuples(), start=1):
    region_name_fig = f"Cluster DOWN rank {i}: {gene.name}"
    region_name_out = f"down_rank{i:03}"
    plot_all_samples(gene, region_name_fig, region_name_out, replicate_sample_ids, "all_replicate_samples")
    plot_all_samples(gene, region_name_fig, region_name_out, merged_sample_ids, "all_merged_samples")
