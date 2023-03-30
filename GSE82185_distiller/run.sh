#!/usr/bin/env bash
#PJM -g ga17
#PJM -L node=1
#PJM -L rscgrp=regular
#PJM -S
#PJM -m e,r
#PJM --restart

. /work/ga17/share/tools/setup_env.sh
nextflow run -resume distiller.nf -params-file project.yml
