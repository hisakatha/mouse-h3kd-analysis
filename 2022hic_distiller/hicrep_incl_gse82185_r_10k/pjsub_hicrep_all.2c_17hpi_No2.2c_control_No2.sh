#!/usr/bin/env bash
#PJM -g ga17
#PJM -L node=1
#PJM -L rscgrp=short
#PJM -L elapse=6:00:00
#PJM -S
#PJM -m e,r
#PJM --restart
/work/00/ga17/share/tools/R-install/bin/Rscript hicrep.single.R 2c_17hpi_No2 2c_control_No2
