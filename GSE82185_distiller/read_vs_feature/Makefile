SHELL := bash
.DELETE_ON_ERROR:
.SECONDEXPANSION:

all: $$(target)

ZCAT := zcat
BEDTOOLS := /work/ga17/share/tools/bedtools2/bin/bedtools
FAI := /work/ga17/share/references/ensembl_mm10/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.fai
MERGE := ./merge_hic_pairs.sh
GSE82185_late2C.pairs.gz:
	mkdir tmpdir_$@ && cd tmpdir_$@ && \
		$(MERGE) GSE82185_late2C ../../results/pairs_library/GSE82185_late2C_No*.mm10.nodups.pairs.gz && \
		mv $@* ../ && cd ../ && rmdir tmpdir_$@

GSE82185_late2C.bedpe.gz: GSE82185_late2C.pairs.gz
%.bedpe.gz: %.pairs.gz
	zcat $< | grep -v "^#" | awk -v OFS="\t" '{print $$2,$$3-1,$$3,$$4,$$5-1,$$5,NR,".",$$6,$$7}' | gzip -c > $@

read_bedpe += GSE82185_late2C.ext500.bedpe.gz
read_bedpe += GSE82185_late2C.ext1000.bedpe.gz
read_bedpe += GSE82185_late2C.ext2500.bedpe.gz

$(read_bedpe): %.bedpe.gz: $$(shell echo $$* | sed -E 's/\.ext[0-9]+//').bedpe.gz
	EXT=$$(echo $@ | sed -E 's/.*ext([0-9]+)\.bedpe\.gz/\1/'); \
	paste <($(ZCAT) $< | cut -f 1-3 | $(BEDTOOLS) slop -b $$EXT -g $(FAI)) \
	      <($(ZCAT) $< | cut -f 4- | $(BEDTOOLS) slop -b $$EXT -g $(FAI)) | gzip -c > $@

target += $(read_bedpe)
