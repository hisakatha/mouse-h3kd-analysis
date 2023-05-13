#!/usr/bin/env bash
# Setup this code repository
set -u
PREFIX="../."
rsync -avRm\
    -f "+ /*/"\
    -f "+ genome_view2/"\
    -f "+ diff_cpm/"\
    -f "- .*"\
    -f "+ *.sh"\
    -f "+ *.py"\
    -f "+ *.R"\
    -f "+ *.xml"\
    -f "+ *.jbrowse"\
    -f "+ Makefile"\
    -f "+ *.makefile"\
    -f "+ README*"\
    -f "+ VERSION*"\
    -f "- *"\
    $@ "$PREFIX/mouse_rnaseq_h3kd" .

rsync -avRm\
    -f "+ /*/"\
    -f "+ tad_strength_back1_w_balancing/"\
    -f "+ tad_strength_back2_w_balancing_wo_diagonal/"\
    -f "- .*"\
    -f "+ *.sh"\
    -f "+ *.py"\
    -f "+ *.R"\
    -f "+ *.xml"\
    -f "+ *.jbrowse"\
    -f "+ Makefile"\
    -f "+ *.makefile"\
    -f "+ README*"\
    -f "+ VERSION*"\
    -f "- *"\
    $@ "../mouse_hic_h3kd/./2022hic_distiller" .
