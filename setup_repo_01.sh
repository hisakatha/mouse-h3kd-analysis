#!/usr/bin/env bash
# Setup this code repository
# Sync from distiller-nf projects
set -u
PREFIX="/work/ga17/share/."
rsync -avRm\
    -f "+ /*/"\
    -f "+ configs/"\
    -f "+ hicrep/"\
    -f "+ read_vs_feature/"\
    -f "- .*"\
    -f "+ *.sh"\
    -f "+ *.py"\
    -f "+ *.nf"\
    -f "+ *.config"\
    -f "+ *.yml"\
    -f "+ Makefile"\
    -f "+ README*"\
    -f "+ VERSION*"\
    -f "- *"\
    $@ "$PREFIX/GSE82185_distiller" .

rsync -avRm\
    -f "+ /*/"\
    -f "+ configs/"\
    -f "+ hicrep/"\
    -f "+ hicrep_incl_gse82185_r_10k/"\
    -f "+ read_vs_feature/"\
    -f "- .*"\
    -f "+ *.sh"\
    -f "+ *.py"\
    -f "+ *.nf"\
    -f "+ *.config"\
    -f "+ *.yml"\
    -f "+ Makefile"\
    -f "+ README*"\
    -f "+ VERSION*"\
    -f "- *"\
    $@ "$PREFIX/2022hic_distiller" .
