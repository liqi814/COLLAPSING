library(tidyverse)
library(tidyr)
library(vcfR)
library(broom)
library(openxlsx)
library(data.table)

readVcf <- function(vcf_file) {
    vcf <- read.vcfR(vcf_file)
    v <- vcfR2tidy(vcf, format_fields = "GT")
    csqd <- v$meta %>% filter(ID == "CSQ") %>% select(Description)
    csqd <- unlist(str_split(gsub("Consequence annotations from Ensembl VEP. Format: ", "", csqd$Description), "\\|"))
    annd <- v$meta %>% filter(ID == "ANN") %>% select(Description)
    annd <- unlist(str_split(gsub("'", "", gsub(" ","" ,gsub("Functional annotations: ", "", annd$Description))), "\\|"))


    CSQ <- v$fix %>%
        separate_rows(CSQ, sep = ",") %>%
        separate(CSQ, csqd, sep = "\\|", convert = TRUE) %>%
                filter(str_detect(Feature, "ENST")) %>%
        unite("VariantID", c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>%
        unite("key", ChromKey, POS, sep = "-", remove = FALSE) %>%
                select(-c(ANN))


    ANN <-  v$fix %>%
        unite("VariantID", c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>%
        select("VariantID", "ANN") %>%
        separate_rows(ANN, sep = ",") %>%
        separate(ANN, annd, sep = "\\|", convert = TRUE) %>%
        filter(str_detect(Feature_ID, "ENST")) %>%
        mutate(Feature_ID = str_extract(as.character(Feature_ID), "ENST[0-9]+"))

    vf <- merge(ANN, CSQ, by.x=c("VariantID", "Feature_ID"), by.y=c("VariantID","Feature"), all.y =TRUE, all.x=FALSE)
    # right join


    vf <- vf %>% mutate(`Variant Type` = case_when(VARIANT_CLASS == "SNV" ~ "snv",
                                          VARIANT_CLASS %in% c("deletion", "insertion") ~ "indel",
                                           TRUE ~ "other")) %>%
        mutate(`Polyphen_char` = case_when(str_detect(PolyPhen, "benign") ~ "benign",
                                           str_detect(PolyPhen, "possibly_damaging") ~ "possibly",
                                           str_detect(PolyPhen, "probably_damaging")  ~ "probably",
                                                   TRUE ~ "unknown")) %>%
        mutate(REVEL_score = as.numeric(REVEL_score),
               PrimateAI_score = as.numeric(PrimateAI_score),
               gnomADe_rf_tp_probability = as.numeric(gnomADe_rf_tp_probability),
               gnomADg_rf_tp_probability = as.numeric(gnomADg_rf_tp_probability),
               ExAC_VQSLOD = as.numeric(ExAC_VQSLOD)) %>%
        mutate(LOO_AF = (as.numeric(AC)-2)/(as.numeric(AN)-2)) %>%
        mutate_at(vars(contains("_AF_") & !contains("MAX")), as.numeric)


    vfp <- vf %>% filter(FILTER == "PASS" & SYMBOL_SOURCE == "HGNC")

    gt <- v$gt %>% unite("key", ChromKey, POS, sep = "-", remove = FALSE)


    return(list(vfp, gt))
}


filterVcf <- function(vcf_data, effect, effect_impact,
                      exac_vqslod_snv, exac_vqslod_indel, 
                      gnomAD_rf_snv, gnomAD_rf_indel,
                      pop, ext_af, loo_af) {

    vcf_fix <- vcf_data[[1]]
    gt <- vcf_data[[2]]
    
    pop2 <- pop[!pop == "asj"]
    
    if (!is.null(effect)) {
        vcf_filt <- vcf_fix %>% 
		filter(str_detect(Consequence, paste(effect, collapse = "|")))

    } else {
        vcf_filt <- vcf_fix
    }
    
    vcf_filt <- vcf_filt %>% 
        group_by(key) %>%
        filter(n_distinct(VariantID) < 2) %>%
        ungroup() %>%
        filter(BIOTYPE == "protein_coding") %>%
        filter(IMPACT %in% effect_impact) %>%
        filter(is.na(`Annotation_Impact`) | `Annotation_Impact` == "" | Annotation_Impact %in% effect_impact) %>%
        filter(CCDS != "" & !(is.na(CCDS))) %>% 
        filter(LOO_AF <= loo_af) %>% 
        filter(is.na(`ExAC_VQSLOD`) | `ExAC_VQSLOD` == "" |
                   `Variant Type` == "snv" & `ExAC_VQSLOD` >= exac_vqslod_snv | 
                   `Variant Type` == "indel" & `ExAC_VQSLOD` >= exac_vqslod_indel) %>%
        filter(is.na(`gnomADe_rf_tp_probability`) | `gnomADe_rf_tp_probability` == "" |
                   `Variant Type` == "snv" & `gnomADe_rf_tp_probability` >= gnomAD_rf_snv | 
                   `Variant Type` == "indel" & `gnomADe_rf_tp_probability` >= gnomAD_rf_indel) %>%
        filter(is.na(`gnomADg_rf_tp_probability`) | `gnomADg_rf_tp_probability` == "" |
                   `Variant Type` == "snv" & `gnomADg_rf_tp_probability` >= gnomAD_rf_snv | 
                   `Variant Type` == "indel" & `gnomADg_rf_tp_probability` >= gnomAD_rf_indel) %>%
        filter(`gnomADg_non_topmed_AF` <= ext_af | is.na(`gnomADg_non_topmed_AF`) | `gnomADg_non_topmed_AF` == "") %>% 
        filter_at(paste0("gnomADe_non_topmed_AF_", pop), all_vars(. <= ext_af | is.na(.) | . == "")) %>% 
        filter_at(paste0("ExAC_AF_", toupper(pop2)), all_vars(. <= ext_af | is.na(.) | . == "")) %>%
        mutate(IMPACT = factor(IMPACT, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))) %>%
        group_by(`VariantID`) %>% 
        arrange(IMPACT) %>%
        distinct(VariantID, .keep_all = TRUE) %>% 
        ungroup()
    
    
    variants <- vcf_filt %>% left_join(gt, by = "key") %>% select(Indiv, everything()) %>% select(-contains(".y"))
    variants <- variants %>% 
        separate(gt_GT_alleles, c("gt_allele1", "gt_allele2"), sep = "/|\\|") %>% 
        mutate(gt_score = case_when(gt_GT %in% c("1/1", "1|1") & (gt_allele1 == ALT & gt_allele2 == ALT) ~ 2,
                                    TRUE ~ 0))
									
    hom_key <- variants %>% group_by(key) %>%
	summarise(no_hom_position = all(gt_score == 0)) %>%
	filter(!no_hom_position) %>% 
	select(key) %>% ungroup()
	
	
    variants <- variants %>% filter(key %in% hom_key$key )
									
    return(variants)
}

collapseVariants <- function(variants, sample_file) {
    print(sample_file)
    ped <- fread(sample_file, header = FALSE, sep = "\t")
    samples <- ped$V2
    is.case <- as.matrix(ped[, 6] == 2)
    rownames(is.case) <- samples

    n.samples <- length(is.case)
    n.cases <- length(which(is.case))
    n.ctrls <- n.samples - n.cases

    if(dim(variants)[1] == 0){
        data <- matrix(0L, nrow = 1, ncol = n.samples)
        rownames(data) <- "" 
        colnames(data) <- samples
	gsum <- data.frame("Gene Name"="",
                           "Qualified Case"=as.double(0),
                           "Unqualified Case"=as.double(n.cases),
                           "Qualified Ctrl"=as.double(0),
                           "Unqualified Ctrl"=as.double(n.ctrls),
                           stringsAsFactors=FALSE)	
    } else {
        coll_matrix <- variants %>% group_by(Indiv, SYMBOL) %>% summarize(sum(gt_score))
        coll_matrix <- coll_matrix %>% pivot_wider(names_from = Indiv, values_from = `sum(gt_score)`)

        data <- coll_matrix %>% column_to_rownames("SYMBOL") %>% as.matrix()

        data[data > 0] <- 1  ###_matrix.txt.gz file
		
        qv.cases <- rowSums(data[,rownames(is.case)[is.case], drop = FALSE])
        qv.ctrls <- rowSums(data[,rownames(is.case)[!is.case], drop = FALSE])

        gsum <- data.frame(qv.cases, n.cases - qv.cases, qv.ctrls, n.ctrls - qv.ctrls, stringsAsFactors = FALSE)
        gsum <- gsum %>% rownames_to_column(var = "Gene Name")}

    colnames(gsum) <- c("Gene Name", "Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl")
    
    return(list(data, gsum))
}

getCollapseData <- function(model, cluster, chr) {
    sample_file <- paste0("/nfs/projects/dbGap/COPD_sample_remove/flashpca/flashPCA_lclustering_res_0_3_cluster_", cluster, "_sample.txt")
    #sample_file <- paste0("/nfs/projects/dbGap/collapsetest/dbGap_greedy_lclustering_res_0_2_cluster_", cluster, "_sample.txt")
    vcf_file <- paste0("/nfs/projects/dbGap/COPD_sample_remove/chr_vcf_cluster/chr_", chr, "/dbGap_all_chr", chr, "_coding_phase2_VEP_cluster_", cluster, "_nomono.vcf.gz")

    results_path <- paste0("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/LClust_res_0_3_cluster_", cluster, "_All_FlashColl_07")
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
    
    if (model == "recessiveAutosomal") {
        effect <- NULL
        effect_impact <- c("HIGH", "MODERATE")
        exac_vqslod_snv = -2.632
        exac_vqslod_indel = 1.262
        gnomAD_rf_snv = 0.01
        gnomAD_rf_indel = 0.02
        pop <- c("afr", "amr", "asj", "eas", "sas", "fin", "nfe")
        ext_af = 0.01
        loo_af = 0.01
    } else {
	stop("ERROR: This script is for recessive model. Please enter model recessiveAutosomal!")
	}
    
    vcf_data <- readVcf(vcf_file)
    variants <- filterVcf(vcf_data, effect, effect_impact,
	                  exac_vqslod_snv, exac_vqslod_indel, gnomAD_rf_snv,
                          gnomAD_rf_indel, pop, ext_af, loo_af)
	
	
    ### add chohort info 
    dbGap_control_IDs <- read.csv(file="/nfs/projects/dbGap/reference/dbGap_control_IDs.csv",header=TRUE) %>%
	    mutate(Sample_ID = as.character(Sample_ID))
    
    qv_variants <- variants %>% 
        left_join(dbGap_control_IDs, by = c("Indiv" = "Sample_ID")) %>% select(cohort, everything()) %>%
        filter(gt_score > 0) %>%
        mutate(cohort = case_when(is.na(cohort) ~ "IPF", TRUE ~ as.character(cohort)))
   

    collap_data <- collapseVariants(variants, sample_file)
	
    data <- collap_data[[1]]
    gsum <- collap_data[[2]]

    if(dim(qv_variants)[1] == 0){
        print(paste0(Sys.Date(),": For model ",  model,  ", No qualified variant found on chr", chr, " in cluster ", cluster))
    } else {
        saveRDS(qv_variants, paste0(model_path, "/", Sys.Date(), "_", model, "_chr_", chr, "_genotypes.RDS"))
        saveRDS(data, paste0(model_path, "/", Sys.Date(), "_", model, "_chr_", chr, "_matrix_txt.RDS"))}
	
    return(gsum)
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

getCMHres <- function(model, chr, nclust) {
    print(model)
    gsum <- lapply(nclust, function(x) getCollapseData(model, x, chr))

    names(gsum) <- as.character(nclust)
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    gene.names.r <- gene.names.r[gene.names.r != ""]
    
    gsum <- lapply(gsum, function(x) x %>% full_join(as_tibble(gene.names.r), by = c("Gene Name" = "value")) %>% 
                       mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
                              `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
                              `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
                              `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
                       )%>% 
                       filter(`Gene Name` != ""))    
    
    cmhs <- lapply(gene.names.r, function(i) useCMH(gsum, i))
    
    cmhs.tidy <- unnest(tibble(cmhs), cols = c(cmhs)) %>% 
        mutate(gene = gene.names.r) %>% arrange(p.value) %>% 
        select(-c(method, alternative, statistic)) %>% 
        select(gene, p.value, everything())
	
    return(cmhs.tidy)
}

##setting##
setwd("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/")
nclust <- c(0:6)
resolution <- "0_3"

args = commandArgs(trailingOnly=TRUE)
model = args[1]
chr = as.numeric(args[2])

###main##
cmhs.tidy <- getCMHres(model, chr, nclust)
rds_path <- paste0("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/", model, "_RDS_path/")
saveRDS(cmhs.tidy, paste0(rds_path, Sys.Date(), "_", model, "_chr", chr, "_res_", resolution,"_cluster_", paste(nclust, collapse = "_"), "_exact_CMH.RDS"))

##main##
#args = commandArgs(trailingOnly=TRUE)
#cluster = args[1]
#getCollapseData(model,cluster , chr)
#lapply(nclust, function(x) getCollapseData(model, x, chr))
