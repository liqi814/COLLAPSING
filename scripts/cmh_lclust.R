library(tidyverse)
library(tidyr)
library(broom)
library(openxlsx)
library(data.table)


## matrix.txt.gz 
combineMatrixRDS <- function(cluster, model) {
    results_path <- paste0("LClust_res_0_3_cluster_", cluster, "_All_FlashColl_07", "/", model)
	
	AUTOSOMES <- c(1:22)
    matrix_path <- lapply(AUTOSOMES, function(x) list.files(results_path, paste0("*chr_", x, "_matrix_txt.RDS$"), full.names = TRUE)[1])
    matrix_file <- lapply(matrix_path, function(x) readRDS(x))
    matrix_files <- do.call("rbind", matrix_file)

	gz1 <- gzfile(paste0(results_path, "/", Sys.Date(), "_", model,"_matrix.csv.gz"), "w")
    write.csv(matrix_files, gz1, quote=FALSE, row.names=TRUE)
	
	close(gz1) 
}

###setting##
setwd("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/")
nclust <- c(0:6)
#model <- "Synonymous"
model <- "RareEnsemble"
#model <- "PTV"
resolution <- "0_3"


dir <- "All"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

CMH_dir <- paste0("CMH_", gsub("_", "",  dir), "/")
CMH_dir <- gsub("_/", "/", CMH_dir)

### combine cmh.tidy 
rds_path <- paste0(model, "_RDS_path/")
cmh_path <- lapply(c(1:22), function(chr) list.files(rds_path, paste0("*chr", chr, "_res_", resolution,"_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS$"), full.names = TRUE)[1])
cmh.list <- lapply(cmh_path, function(x) readRDS(x))
cmh.all <- do.call("rbind", cmh.list)  %>% arrange(p.value)

write.xlsx(cmh.all, paste0(CMH_dir, "All_", model,"_CMH_exact_summary_lclust_res_", resolution,".xlsx"), row.names=FALSE)


### combine matrix files and prepare for QQ plot
lapply(nclust,function(cluster) combineMatrixRDS(cluster, model))