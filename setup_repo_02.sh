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
