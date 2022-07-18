#! /bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -pe threaded 24
source activate ATAV
Rscript /nfs/projects/dbGap/COPD_sample_remove/flashpca/scripts/lclust_Flash.R 0.3
