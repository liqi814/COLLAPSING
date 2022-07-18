#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-02-01 11:59
# * Last modified : 2021-02-01 11:59
# * Filename      : qsub_cmh_permute_lclust_short.sh
# * Description   :


LOG_DIR="/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/LatinoClustersAll/Log"

for NAME in Synonymous 
do
  OUT=$LOG_DIR/$NAME\_permute_stdout.log
  ERR=$LOG_DIR/$NAME\_permute_stderr.log

  QSUB_CMD="qsub -cwd -V -S /bin/bash -N $NAME  -pe threaded 24 -o $OUT -e $ERR"
  ATAV_CMD="source activate ATAV; Rscript /nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/scripts/cmh_permute_lclust_short.R $NAME"

  echo $ATAV_CMD | $QSUB_CMD
done
