#!/bin/bash

# * Author        : Qi Li
# * Email         : ql2387@cumc.columbia.edu
# * Create time   : 2021-01-11 13:43
# * Last modified : 2021-01-11 13:43
# * Filename      : split_cluster_no_mono.sh
# * Description   :

#for CHROM in 3 4 $(seq 7 22)
#for CHROM in 2 5
#for CHROM in 1 $(seq 3 9) 11 $(seq 13 22)
for CHROM in $(seq 1 22)
do

mkdir -p "/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/chr_"$CHROM
LOGS="/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/logs"

for CLUSTER in {0..6}
do

BASH_SCRIPTS="/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/scripts/chr_"$CHROM"_cluster_"$CLUSTER"_no_mono_wLOFTEE_IMP.sh"

if [ -f "$BASH_SCRIPTS" ]; then
    rm -f $BASH_SCRIPTS
fi

echo "source /nfs/goldstein/software/bcftools-1.9-x86_64/bcftools-1.9-ENV.sh " >> $BASH_SCRIPTS

echo SAMPLE_list="/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/sample_list_1c/flashPCA_lclustering_res_0_3_cluster_"$CLUSTER"_sample_1c.txt" >> $BASH_SCRIPTS
echo CLUSTER_OUTPUT="/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/chr_"$CHROM"/dbGap_all_chr"$CHROM"_coding_phase2_VEP_LOFTEE_IMP_cluster_"$CLUSTER".vcf.gz" >> $BASH_SCRIPTS
echo bcftools view -S \$SAMPLE_list -Oz "/nfs/projects/dbGap/dbGap_annotated/dbGap_all_chr"$CHROM"_coding_phase2_VEP_LOFTEE_IMP*.vcf.gz" \> \$CLUSTER_OUTPUT  >> $BASH_SCRIPTS
echo NO_MONO_CLUSTER_OUTPUT="/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/chr_"$CHROM"/dbGap_all_chr"$CHROM"_coding_phase2_VEP_LOFTEE_IMP_cluster_"$CLUSTER"_nomono.vcf.gz" >> $BASH_SCRIPTS
echo bcftools index -t -f \$CLUSTER_OUTPUT >> $BASH_SCRIPTS
echo bcftools view -c 1:'nonmajor' -Oz \$CLUSTER_OUTPUT \> \$NO_MONO_CLUSTER_OUTPUT  >> $BASH_SCRIPTS

sh /home/ql2387/scripts/0.qsub_jobs_directory.sh  $(basename $BASH_SCRIPTS .sh) $BASH_SCRIPTS  $LOGS

done
done 
