#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-01-07 12:10
# * Last modified : 2021-01-27 13:03
# * Filename      : kinship_and_pruning.sh
# * Description   :

##kinship analysis
#/nfs/goldstein/software/king_relatedness/king \
#-b /nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_all_kin_snps_hg19_exact_nodup.bed \
#--kinship --related --degree 3 \
#--prefix /nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_noIGMdup_noCOPD

##kinship pruning
/nfs/goldstein/software/python2.7.7/bin/python \
/nfs/goldstein/software/atav_home/lib/run_kinship.py \
/nfs/projects/dbGap/COPD_sample_remove/KING/dbGap_noIGM_noCOPD_sample.txt \
/nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_noIGM_noCOPD.kin0 \
/nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_noIGM_noCOPD.kin \
--relatedness_threshold 0.0884 --seed 42 \
--output /nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_noIGM_noCOPD_pruned_sample.txt
