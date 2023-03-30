#!/usr/bin/env bash
set -ue
# First argument is an output prefix.
# The rest arguments are input .pairs.gz files.
# OUTPUT_PREFIX.pairs.gz will be created.
PAIRIX_ROOT="/work/ga17/share/tools/pairix-0.3.7"
export PATH="$PAIRIX_ROOT/bin:$PATH"
"$PAIRIX_ROOT/util/merge-pairs.sh" "$@"
