library(tidyverse)
library(tidyr)
library(broom)
library(openxlsx)
library(data.table)


###setting##
#model <- "recessiveAutosomal"
#model <- "PTVwLOFTEE"
setwd("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/")
model <- "RareEnsemble"
#model <- "PTV"
#model <- "dominateControl"
resolution <- "0_3"
args = commandArgs(trailingOnly=TRUE)
chr <- as.numeric(args[1])
GENE_NAME <- args[2]
OUT_DIR <- paste0(GENE_NAME, "_variants/")

if(!dir.exists(OUT_DIR)) {
        dir.create(OUT_DIR)
    } else {
        print((paste0("Directory ", OUT_DIR, " already exists")))
    }


nclust <- c(3)
#genotypes
#--list-geno-var
genotypes_files <- lapply(nclust, function(x) list.files(paste0("LClust_res_0_3_cluster_", x, "_All_FlashColl_07", "/", model), paste0("*chr_",chr ,"_genotypes.RDS$"), full.names = TRUE)[1])
genotypes <- lapply(genotypes_files, function(x) if (!is.na(x)) {readRDS(x)})
#genotypes <- lapply(genotypes_files, function(x) readRDS(x))
genotype_all <- rbindlist(genotypes)

QV <- genotype_all  %>% filter(SYMBOL == GENE_NAME) %>% filter(gt_score >0)
#write.xlsx(QV, paste0(OUT_DIR, Sys.Date(),"_", model, "_", GENE_NAME, "_filter_variants_nlcluster.xlsx"), row.names=FALSE)
write.csv(QV, paste0(OUT_DIR, Sys.Date(),"_", model, "_", GENE_NAME, "_filter_variants_nlcluster.csv"), row.names=FALSE)

##whole chromosome
#genotypes_files <- lapply(1:22, function(CHR) lapply(nclust, function(x) list.files(paste0("LClust_res_0_3_cluster_", x, "_All_FlashColl_07", "/", model), paste0("*chr_",CHR ,"_genotypes.RDS$"), full.names = TRUE)[1]))
#genotypes_files <- unlist(genotypes_files, recursive = FALSE)

