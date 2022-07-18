library(tidyverse)
library(tidyr)
library(broom)
library(openxlsx)
library(data.table)

setwd("/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/")

getCollapseData <- function(model, cluster, res) {
    sample_file <- paste0("flashpcaAll/flashPCA_lclustering_res_", res, "_cluster_", cluster, "_sample.txt")
    results_path <- paste0("CollapsingAll/LClust_res_", res, "_cluster_", cluster, "_All_FlashColl_07")
	
    if(!dir.exists(results_path)) {
        dir.create(results_path)
    } else {
        print((paste0("Directory ", results_path, " already exists")))
    }

    model_path <- paste0(results_path, "/", model)
    if(!dir.exists(model_path)) {
        dir.create(model_path)
    } else {
        print((paste0("Directory ", model_path, " already exists")))
    }

    ped <- fread(sample_file, header = FALSE, sep = "\t")
    samples <- ped$V2
    is.case <- as.matrix(ped[, 6] == 2)
    rownames(is.case) <- samples

    n.samples <- length(is.case)
    n.cases <- length(which(is.case))
    n.ctrls <- n.samples - n.cases
	
    genotype_file <- list.files(model_path, "_genotypes.csv.gz$", full.names = TRUE)
	
    variants <- fread(genotype_file)
    variants$`Gene Name` <- variants$`Gene Name` %>%  str_replace_all("'", "")
    variants$`Sample Name` <- variants$`Sample Name` %>%  str_replace_all("'", "")
    variants$GT <- variants$GT %>%  str_replace_all("'", "")

    variants <- variants %>% 
    mutate(gt_score = case_when(GT == "het" ~ 1,
                                GT == "hom" ~ 2,
                                TRUE ~ 0))								

    coll_matrix <- variants %>% group_by(`Sample Name`, `Gene Name`) %>% summarize(sum(gt_score))
    coll_matrix <- coll_matrix %>% pivot_wider(names_from = `Sample Name`, values_from = `sum(gt_score)`, values_fill = 0)

    data <- coll_matrix %>% column_to_rownames("Gene Name") %>% as.matrix()

    if (dim(data)[2] < n.samples) {
        data <- cbind(data,0)
        colnames(data)[dim(data)[2]] <- setdiff(rownames(is.case), colnames(data))}

    data[data > 0] <- 1  ###_matrix.txt.gz file

    gz1 <- gzfile(paste0(model_path, "/", Sys.Date(), "_", model,"_matrix.csv.gz"), "w")
    write.csv(data, gz1, quote=FALSE, row.names=TRUE)
    close(gz1)
	
    qv.cases <- rowSums(data[,rownames(is.case)[is.case], drop = FALSE])
    qv.ctrls <- rowSums(data[,rownames(is.case)[!is.case], drop = FALSE])

    gsum <- data.frame(qv.cases, n.cases - qv.cases, qv.ctrls, n.ctrls - qv.ctrls, stringsAsFactors = FALSE)
    gsum <- gsum %>% rownames_to_column(var = "Gene Name")
    
    colnames(gsum) <- c("Gene Name", "Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl")

    saveRDS(gsum, paste0(model_path, "/", Sys.Date(), "_", model, "_exact_CMH.RDS"))	
}

#main 
model <- "dominantRareEnsemble2New2"
nclust <- c(0:10)
resolution <- "0_25"

lapply(nclust, function (x) getCollapseData(model, x, resolution))
