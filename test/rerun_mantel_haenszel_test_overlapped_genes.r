library(tidyverse)
library(tidyr)
library(vcfR)
library(broom)
library(openxlsx)
library(data.table)

read_excel_allsheets <- function(filename, tibble = FALSE) {
    # I prefer straight data.frames
    # but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    x
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

getCMHres <- function(model, IGM_nclust, IGM_out, dbGAP_nclust) {
    print(model)
	
	model_outputs <- readData(model, IGM_out, dbGAP_nclust)
	
	gsum <- getgsumData(IGM_nclust, dbGAP_nclust, model_outputs)
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    
    #gsum <- lapply(gsum, function(x) x %>% full_join(as_tibble(gene.names.r), by = c("Gene Name" = "value")) %>% 
    #                   mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
    #                          `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
    #                          `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
    #                          `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
    #                   ))    
    
    cmhs <- lapply(gene.names.r, function(i) useCMH(gsum, i))
    
    cmhs.tidy <- unnest(tibble(cmhs), cols = c(cmhs)) %>% 
        mutate(gene = gene.names.r) %>% arrange(p.value) %>% 
        select(-c(method, alternative, statistic)) %>% 
        select(gene, p.value, everything())   

    write.xlsx(cmhs.tidy, paste0("All_IGM_dbGAP_",model, "_cluster_exact_CMH.xlsx"), rownames=FALSE)
	return(cmhs.tidy)
}


getgsumData <- function(IGM_nclust, dbGAP_nclust, model_outputs) {
	gsum_colnames <- c("Gene Name","Qualified Case","Unqualified Case","Qualified Ctrl","Unqualified Ctrl")
	gsum <- list()
	index=0
	
	##keep only overlapped gene
	IGM_dbGAP_model_out <- merge(model_outputs[[1]], model_outputs[[2]], by="gene")
	
	for (clst in IGM_nclust)
	{
		index=index+1
		gsum[[index]] <- IGM_dbGAP_model_out %>% 
			select(gene, paste(clst,".Qualified Case", sep=""), paste(clst,".Unqualified Case", sep=""),paste(clst,".Qualified Ctrl", sep=""),paste(clst,".Unqualified Ctrl", sep=""))
		colnames(gsum[[index]]) <- gsum_colnames
	}
	for (clst in dbGAP_nclust)
	{
		index=index+1
		gsum[[index]] <- IGM_dbGAP_model_out %>% 
			select(gene, paste(clst,".Qualified.Case", sep=""), paste(clst,".Unqualified.Case", sep=""),paste(clst,".Qualified.Ctrl", sep=""),paste(clst,".Unqualified.Ctrl", sep=""))
		colnames(gsum[[index]]) <- gsum_colnames
	}
	
	names(gsum) <- as.character(c(paste("IGM" ,IGM_nclust,sep="."), paste("dbGAP",dbGAP_nclust,sep=".")))
	return(gsum)
}


readData <- function(model, IGM_out, dbGAP_nclust){
    resolution <- "0_3"
	CMH_path <- ("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/CMH_All/")
	
	IGM_model_out = IGM_out[[grepl(model, names(IGM_out))]]
	IGM_model_out$gene <- gsub("'(.*)'", "\\1", IGM_model_out$gene)
    
	#dbGAP_model_out = read.xlsx(paste0(CMH_path, "All_", model, "_res_", resolution,"_cluster_", paste(dbGAP_nclust, collapse = "_"), "_exact_summary_CMH.xlsx"))
        dbGAP_model_out = read.xlsx(paste0(CMH_path, "All_", model, "_CMH_exact_summary_lclust_res_",resolution ,".xlsx"))
	
	return(list(IGM_model_out, dbGAP_model_out))
}

model <- "recessiveAutosomal"
dbGAP_nclust <- c(0:6)
IGM_nclust <- c(0:10)

###read IGM collapsing analyses result and pre-process the dataframe
IGM_out <- read_excel_allsheets("/nfs/projects/dbGap/entire_cohort_IGM_CA/Gundula_RE2_Output_2.12.21/All_CMH_exact_summary_lclust_res_0_25_min_sample_5_recessive.xlsx")
#IGM_out <- lapply(IGM_out, function(x) gsub("'(.*)'", "\\1", x$gene))

cmhtidy <- getCMHres(model, IGM_nclust, IGM_out, dbGAP_nclust)

#QV_dbGAP_RareEnsemble_out <- merge(QV_RareEnsemble_out, dbGAP_RareEnsemble_out, by="gene")
#QV_dbGAP_RareEnsemble_out %>% select("Qualified Case", "Qualified.Case") %>% rowSums(na.rm=TRUE) -> QV_dbGAP_RareEnsemble_out$Toal_Qualified_Case
#QV_dbGAP_RareEnsemble_out %>% select("Unqualified Case", "Unqualified.Case") %>% rowSums(na.rm=TRUE) -> QV_dbGAP_RareEnsemble_out$Toal_Unqualified_Case
#QV_dbGAP_RareEnsemble_out %>% select("Qualified Ctrl", "Qualified.Ctrl") %>% rowSums(na.rm=TRUE) -> QV_dbGAP_RareEnsemble_out$Toal_Qualified_Ctrl
#QV_dbGAP_RareEnsemble_out %>% select("Unqualified Ctrl", "Unqualified.Ctrl") %>% rowSums(na.rm=TRUE) -> QV_dbGAP_RareEnsemble_out$Toal_Unqualified_Ctrl

#QV_PTV_out <- IGM_out$PTV
#QV_RareEnsemble_out <- IGM_out$RareEnsemble
#QV_Synonymous_out <- IGM_out$Synonymous



