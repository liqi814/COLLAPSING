#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2022-01-18 14:22
# * Last modified : 2022-01-18 14:22
# * Filename      : list_geno_list.sh
# * Description   : 

GENENAME=$1

EuroIGMnclust="0 1 2 10"
NonEuroIGMnclust="3 4 5 6 7 8 9"
EurodbGAPnclust="0 1 4 6"
NonEurodbGAPnclust="2 3 5"

zcat /nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/CollapsingAll/LClust_res_0_25_cluster_0_All_FlashColl_07/dominantRareEnsemble2New2/*dominantRareEnsemble2New2_genotypes.csv.gz  | head -n 1 > /nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/genotypes_file/${GENENAME}_igm_genotypes_CaseAndCtrl.csv
for CLUSTER in `echo $EuroIGMnclust`
do
zcat  /nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/CollapsingAll/LClust_res_0_25_cluster_${CLUSTER}_All_FlashColl_07/dominantRareEnsemble2New2/*dominantRareEnsemble2New2_genotypes.csv.gz | grep -E \'$GENENAME\' >> /nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/genotypes_file/${GENENAME}_igm_genotypes_CaseAndCtrl.csv

done
