library(tidyverse)
library(tidyr)
library(broom)
library(openxlsx)

getCollapseData <- function(cmhTable, res, cluster){
	paste0(cluster, ".", c("Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl"))
	gsum <- cmhTable %>% select(c("gene", paste0(cluster, ".", c("Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl"))))
	colnames(gsum) <- c("Gene Name", "Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl")
	
	results_path <- paste0("CollapsingAll/LClust_res_", res, "_cluster_", cluster, "_All_FlashColl_07")
	model_path <- paste0(results_path, "/", model)
    if(!dir.exists(model_path)) {
        dir.create(model_path)
    } else {
        print((paste0("Directory ", model_path, " already exists")))
    }

        matrix.file <- list.files(model_path, pattern="matrix.csv.gz$", full.names = TRUE)
        data <- read.csv(matrix.file, sep = ",", row.names=1) %>% as.matrix()
        gsum <- gsum %>% filter(`Gene Name` %in% rownames(data))

	saveRDS(gsum, paste0(model_path, "/", Sys.Date(), "_", model, "_exact_CMH.RDS"))
}

getdbGAPdata <- function(model, nclust, res) {
	model <- gsub("dominant", "", model)
    rds_path <- paste0("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/", model, "_RDS_path/")
	cmh_path <- lapply(c(1:22), function(chr) list.files(rds_path, paste0("*chr", chr, "_res_", res,"_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS$"), full.names = TRUE)[1])
	cmh.list <- lapply(cmh_path, function(x) readRDS(x))
	cmh.all <- do.call("rbind", cmh.list)
	
	lapply(nclust, function(x) getCollapseData(cmh.all, res, x))
}

useCMH <- function(gsum, gene) {
    print(gene)
    conf <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>%
                       select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`))
    confa <- simplify2array(lapply(conf, function(x) array(unlist(x), dim = c(2,2))))

    cmh <- mantelhaen.test(confa, exact = TRUE)

    confp <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>%
                        select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`))

    confb <- do.call(rbind, conf)
    confb <- bind_rows(colSums(confb)) %>%
        mutate(`%QV+ Case` = `Qualified Case`/(`Qualified Case` + `Unqualified Case`)*100,
               `%QV+ Ctrl` = `Qualified Ctrl`/(`Qualified Ctrl` + `Unqualified Ctrl`)*100)

    tcmh <- as_tibble(c(tidy(cmh), confb, unlist(confp)))

    return(tcmh)
}

getCMHres <- function(model, IGMnclust, dbGAPnclust, IGM_resolution, dbGAP_resolution) {
    print(model)
	
	gsumIGMpath <- lapply(IGMnclust, function(cluster) list.files(paste0("CollapsingAll/LClust_res_", IGM_resolution, "_cluster_", cluster, "_All_FlashColl_07"), "_exact_CMH.RDS$", full.names = TRUE, recursive = TRUE)[1])
	gsumdbGAPpath <- lapply(dbGAPnclust, function(cluster) list.files(paste0("CollapsingAll/LClust_res_", dbGAP_resolution, "_cluster_", cluster, "_All_FlashColl_07"), "_exact_CMH.RDS$", full.names = TRUE, recursive = TRUE)[1])
	gsum_path <- c(gsumIGMpath, gsumdbGAPpath)
	
    gsum <- lapply(gsum_path, function(x) readRDS(x))
    names(gsum) <- as.character(c(paste("IGM" ,IGMnclust, sep="."), paste("dbGAP",dbGAPnclust, sep=".")))
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    gene.names.r <- gene.names.r[gene.names.r != ""]

    gsum <- lapply(gsum, function(x) x %>% full_join(as_tibble(gene.names.r), by = c("Gene Name" = "value")) %>%
                       mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
                              `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
                              `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
                              `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
                       ) %>% filter(`Gene Name` != ""))

    cmhs <- lapply(gene.names.r, function(i) useCMH(gsum, i))

    cmhs.tidy <- unnest(tibble(cmhs), cols = c(cmhs)) %>%
        mutate(gene = gene.names.r) %>% arrange(p.value) %>%
        select(-c(method, alternative, statistic)) %>%
        select(gene, p.value, everything())

    cmhs.tidy <-  cmhs.tidy %>% filter(gene != "")
	return(cmhs.tidy)
}


LatinoIGMnclust <- c(3, 6)
IGMnclust <-  c(0:10)

dbGAPnclust <-  c(0:6)
LatinodbGAPnclust <- c(3)

IGM_resolution <- "0_25"
dbGAP_resolution <- "0_3"

model <- "dominantRareEnsemble"
CMH_dir <- "CMH_All/"

# adapt to your path
setwd("/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/")

#mian
#getdbGAPdata(model, dbGAPnclust, dbGAP_resolution)

Latino.cmhs.tidy <- getCMHres(model, LatinoIGMnclust, LatinodbGAPnclust, IGM_resolution, dbGAP_resolution)
write.xlsx(Latino.cmhs.tidy, paste0(CMH_dir, "All_Latino_", model,"_CMH_exact_summary_lclust.xlsx"), row.names=FALSE)

