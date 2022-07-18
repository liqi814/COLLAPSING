#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-01-27 14:39
# * Last modified : 2021-01-27 14:39
# * Filename      : lclust_Flash.sh
# * Description   :

LOGS=/nfs/projects/dbGap/COPD_sample_remove/flashpca/logs

for RESOLUTION in $(seq 0.1 0.1 0.4)
do
  BASH_SCRIPTS="/nfs/projects/dbGap/COPD_sample_remove/flashpca/scripts/lclust_Flash_${RESOLUTION}.sh"
  if [ -f "$BASH_SCRIPTS" ]; then
    rm -f $BASH_SCRIPTS
  fi
  
  echo "#! /bin/bash"  >> $BASH_SCRIPTS
  echo "#$ -S /bin/bash" >> $BASH_SCRIPTS
  echo "#$ -cwd" >> $BASH_SCRIPTS
  echo "#$ -j y" >> $BASH_SCRIPTS
  echo "#$ -V" >> $BASH_SCRIPTS
  echo "#$ -pe threaded 24" >> $BASH_SCRIPTS

  echo "source activate ATAV" >> $BASH_SCRIPTS
  echo "Rscript /nfs/projects/dbGap/COPD_sample_remove/flashpca/scripts/lclust_Flash.R ${RESOLUTION}" >> $BASH_SCRIPTS
  sh /home/ql2387/scripts/0.qsub_jobs_directory.sh  $(basename $BASH_SCRIPTS .sh) $BASH_SCRIPTS  $LOGS
done
