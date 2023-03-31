SHELL := bash

.DELETE_ON_ERROR:
.SECONDEXPANSION:

all: $$(target)

fai := ../../mouse_ref/Ensembl_direct/GRCm38/Sequence/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai
read_template := Aligned.sortedByCoord.out.uniq_collated.template.sorted.bed
read_template_tdf := $(read_template).tdf
target += $(read_template_tdf)
$(read_template_tdf): $(read_template)
	~/softwares/IGV_2.12.2/igvtools count $< $@ $(fai)

read_template_cov_bedgraph := $(read_template).cov.bedGraph
read_template_cov_bigwig := $(read_template).cov.bigWig
target += $(read_template_cov_bigwig)
$(read_template_cov_bigwig): $(read_template)
	bedtools genomecov -bg -i $< -g $(fai) > $(read_template_cov_bedgraph)
	~/softwares/ucsc/bedGraphToBigWig $(read_template_cov_bedgraph) $(fai) $@
	$(RM) $(read_template_cov_bedgraph)

read_template_cpm_bedgraph := $(read_template).cpm.bedGraph
read_template_cpm_bigwig := $(read_template).cpm.bigWig
target += $(read_template_cpm_bigwig)
$(read_template_cpm_bigwig): $(read_template)
	cpm_factor=$$(python3 -c "print(1000000 / $$(cat $< | wc -l))");\
	bedtools genomecov -bg -i $< -g $(fai) -scale $$cpm_factor > $(read_template_cpm_bedgraph)
	~/softwares/ucsc/bedGraphToBigWig $(read_template_cpm_bedgraph) $(fai) $@
	$(RM) $(read_template_cpm_bedgraph)

read_template_cpm_bedgraph_all := $(read_template).cpm.all.bedGraph.gz
target += $(read_template_cpm_bedgraph_all)
$(read_template_cpm_bedgraph_all): $(read_template)
	cpm_factor=$$(python3 -c "print(1000000 / $$(cat $< | wc -l))");\
	bedtools genomecov -bga -i $< -g $(fai) -scale $$cpm_factor | bedtools sort | bgzip -c > $@

BEDMAP := ~/softwares/bedops/bin/bedmap
fai_window_1k := ../../mouse_ref/Ensembl_direct/GRCm38/Sequence/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai.sorted.window_1k.bed
read_template_cpm_window_1k_mean := $(read_template).cpm.window_1k_mean.bedGraph.gz
target += $(read_template_cpm_window_1k_mean)
$(read_template_cpm_window_1k_mean): $(read_template_cpm_bedgraph_all)
	$(BEDMAP) --faster --delim "\t" --echo --wmean $(fai_window_1k) <(zcat $< | sed -E 's/([-.0-9]+)$$/.\t\1/') | bgzip -c > $@
