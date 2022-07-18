#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-02-10 10:43
# * Last modified : 2021-02-10 10:43
# * Filename      : qsub_cmh_qq.sh
# * Description   :

LOG_DIR="/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/logs"

for NAME in Synonymous
do
  OUT=$LOG_DIR/$NAME\_qq_stdout.log
  ERR=$LOG_DIR/$NAME\_qq_stderr.log

  QSUB_CMD="qsub -cwd -V -S /bin/bash -N $NAME  -pe threaded 24 -o $OUT -e $ERR  -hold_jid 4947086"
  ATAV_CMD="source activate ATAV; Rscript /nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/scripts/cmh_qq.R $NAME"

  echo $ATAV_CMD | $QSUB_CMD
done 
