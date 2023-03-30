SHELL := bash
.DELETE_ON_ERROR:
.SECONDEXPANSION:
.PHONY: all

all: $$(target)

ZCAT := zcat
BEDTOOLS := /work/ga17/share/tools/bedtools2/bin/bedtools
FAI := /work/ga17/share/references/ensembl_mm10/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai
MERGE := ./merge_hic_pairs.sh
2c_control.pairs.gz:
	mkdir tmpdir_$@ && cd tmpdir_$@ && \
		$(MERGE) 2c_control ../../results/pairs_library/2c_control_No*.mm10.nodups.pairs.gz && \
		mv $@* ../ && cd ../ && rmdir tmpdir_$@
2c_17hpi.pairs.gz:
	mkdir tmpdir_$@ && cd tmpdir_$@ && \
		$(MERGE) 2c_17hpi ../../results/pairs_library/2c_17hpi_No*.mm10.nodups.pairs.gz && \
		mv $@* ../ && cd ../ && rmdir tmpdir_$@
2c_H3-1_and_3-2KD.pairs.gz:
	mkdir tmpdir_$@ && cd tmpdir_$@ && \
		$(MERGE) 2c_H3-1_and_3-2KD ../../results/pairs_library/2c_H3-1_and_3-2KD_No*.mm10.nodups.pairs.gz && \
		mv $@* ../ && cd ../ && rmdir tmpdir_$@

2c_control.bedpe.gz: 2c_control.pairs.gz
bedpe_from_merged_pairs += 2c_control.bedpe.gz
2c_17hpi.bedpe.gz: 2c_17hpi.pairs.gz
bedpe_from_merged_pairs += 2c_17hpi.bedpe.gz
2c_H3-1_and_3-2KD.bedpe.gz: 2c_H3-1_and_3-2KD.pairs.gz
bedpe_from_merged_pairs += 2c_H3-1_and_3-2KD.bedpe.gz

2c_control_No1.bedpe.gz: ../results/pairs_library/2c_control_No1.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_control_No1.bedpe.gz
2c_control_No2.bedpe.gz: ../results/pairs_library/2c_control_No2.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_control_No2.bedpe.gz
2c_17hpi_No1.bedpe.gz: ../results/pairs_library/2c_17hpi_No1.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_17hpi_No1.bedpe.gz
2c_17hpi_No2.bedpe.gz: ../results/pairs_library/2c_17hpi_No2.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_17hpi_No2.bedpe.gz
2c_H3-1_and_3-2KD_No1.bedpe.gz: ../results/pairs_library/2c_H3-1_and_3-2KD_No1.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_H3-1_and_3-2KD_No1.bedpe.gz
2c_H3-1_and_3-2KD_No2.bedpe.gz: ../results/pairs_library/2c_H3-1_and_3-2KD_No2.mm10.nodups.pairs.gz
bedpe_from_single_pairs += 2c_H3-1_and_3-2KD_No2.bedpe.gz

$(bedpe_from_merged_pairs): %.bedpe.gz: %.pairs.gz
	zcat $< | grep -v "^#" | awk -v OFS="\t" '{print $$2,$$3-1,$$3,$$4,$$5-1,$$5,NR,".",$$6,$$7}' | gzip -c > $@
$(bedpe_from_single_pairs): %.bedpe.gz: ../results/pairs_library/%.mm10.nodups.pairs.gz
	zcat $< | grep -v "^#" | awk -v OFS="\t" '{print $$2,$$3-1,$$3,$$4,$$5-1,$$5,NR,".",$$6,$$7}' | gzip -c > $@

#PAIRS2BEDPE := /work/ga17/share/hic_analysis/pairs2bedpe.sh
read_bedpe += 2c_control.ext500.bedpe.gz
read_bedpe += 2c_control.ext1000.bedpe.gz
read_bedpe += 2c_control.ext2500.bedpe.gz

read_bedpe += 2c_17hpi.ext500.bedpe.gz
read_bedpe += 2c_17hpi.ext1000.bedpe.gz
read_bedpe += 2c_17hpi.ext2500.bedpe.gz

read_bedpe += 2c_H3-1_and_3-2KD.ext500.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD.ext1000.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD.ext2500.bedpe.gz

read_bedpe += 2c_control_No1.ext500.bedpe.gz
read_bedpe += 2c_control_No1.ext1000.bedpe.gz
read_bedpe += 2c_control_No1.ext2500.bedpe.gz
read_bedpe += 2c_control_No2.ext500.bedpe.gz
read_bedpe += 2c_control_No2.ext1000.bedpe.gz
read_bedpe += 2c_control_No2.ext2500.bedpe.gz

read_bedpe += 2c_17hpi_No1.ext500.bedpe.gz
read_bedpe += 2c_17hpi_No1.ext1000.bedpe.gz
read_bedpe += 2c_17hpi_No1.ext2500.bedpe.gz
read_bedpe += 2c_17hpi_No2.ext500.bedpe.gz
read_bedpe += 2c_17hpi_No2.ext1000.bedpe.gz
read_bedpe += 2c_17hpi_No2.ext2500.bedpe.gz

read_bedpe += 2c_H3-1_and_3-2KD_No1.ext500.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD_No1.ext1000.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD_No1.ext2500.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD_No2.ext500.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD_No2.ext1000.bedpe.gz
read_bedpe += 2c_H3-1_and_3-2KD_No2.ext2500.bedpe.gz
$(read_bedpe): %.bedpe.gz: $$(shell echo $$* | sed -E 's/\.ext[0-9]+//').bedpe.gz
	EXT=$$(echo $@ | sed -E 's/.*ext([0-9]+)\.bedpe\.gz/\1/'); \
	paste <($(ZCAT) $< | cut -f 1-3 | $(BEDTOOLS) slop -b $$EXT -g $(FAI)) \
	      <($(ZCAT) $< | cut -f 4- | $(BEDTOOLS) slop -b $$EXT -g $(FAI)) | gzip -c > $@

target += $(read_bedpe)

bams := 2c_17hpi.bam 2c_17hpi_No1.bam 2c_17hpi_No2.bam
bams += 2c_control.bam 2c_control_No1.bam 2c_control_No2.bam
bams += 2c_H3-1_and_3-2KD.bam 2c_H3-1_and_3-2KD_No1.bam 2c_H3-1_and_3-2KD_No2.bam
$(bams): %.bam: %.bedpe.gz
	$(BEDTOOLS) bedpetobam -i $< -g $(FAI) > $@
target += $(bams)

sorted_bams := $(bams:%.bam=%.sorted.bam)
$(sorted_bams): %.sorted.bam: %.bam
	samtools sort -o $@ -@ 3 $<

bais := $(sorted_bams:%=%.bai)
$(bais): %.bai: %
	samtools index $<

target += $(bais)
