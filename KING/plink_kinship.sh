#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-01-07 12:06
# * Last modified : 2021-01-07 12:06
# * Filename      : plink_kinship.sh
# * Description   :

plink --bfile /nfs/projects/dbGap/kin_flashpca/dbgap_prep/dbGap_all_kin_snps_hg19_extact --make-bed \
--out /nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_all_kin_snps_hg19_exact_nodup \
--remove /nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/cross_dup_199_plus_COPD.txt 
