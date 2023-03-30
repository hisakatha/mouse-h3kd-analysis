#!/usr/bin/env bash
echo "Do not use this"
echo "coolbox image at docker://nanguage/coolbox is obsolete"
exit

module load singularity
singularity exec docker://nanguage/coolbox ./coolbox_plot.py
