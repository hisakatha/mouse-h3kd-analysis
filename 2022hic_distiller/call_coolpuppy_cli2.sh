#!/usr/bin/env bash
pileup=GSE82185_late2C.ATACloss_ZGA.clpy
coolpup.py GSE82185_late2C.mm10.no_filter.1000.mcool::resolutions/10000 enhancer_tss_loci/ATACloss_ZGA.pairs.ensembl.bedpe --features_format bedpe --nshifts 1 -o $pileup
plotpup.py --input_pups $pileup --no_score --height 2.5 --scale linear --not_symmetric --vmin 0.97 --vmax 1.1 -o $pileup.pdf
