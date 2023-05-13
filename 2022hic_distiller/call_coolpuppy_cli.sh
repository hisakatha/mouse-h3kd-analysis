#!/usr/bin/env bash
pileup=GSE82185_late2C.1000kb_atac_major_zga.clpy
coolpup.py GSE82185_late2C.mm10.no_filter.1000.mcool::resolutions/10000 enhancer_tss_loci/1000kb_atac_major_zga.ensembl.bedpe --features_format bedpe --nshifts 1 -o $pileup
plotpup.py --input_pups $pileup --not_symmetric -o $pileup.pdf
