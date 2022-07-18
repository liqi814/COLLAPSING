# library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

library(data.table)
library(broom)
library(parallel)

# adapt to your data structure
readMatrix <- function(model, dir, pc_path, cluster) {
    sample.file <- paste0(pc_path, "/flashPCA_lclustering_res_", resolution, "_cluster_", cluster, "_sample.txt")
    matrix.file <- list.files(paste0("LClust_res_", resolution, "_cluster_", cluster, "_", dir, "FlashColl_07/", model), "_matrix.csv.gz$", full.names = TRUE)
    
    ped <- fread(sample.file, header = FALSE, sep = "\t")
    samples <- ped$V2
    is.case <- as.matrix(ped[, 6] == 2)
    rownames(is.case) <- samples
    print(matrix.file)
    #    data <- fread(matrix.file, header = TRUE, sep = "\t")
    #data <- read.csv(matrix.file, delim = "\t", row.names=1)
	data <- read.csv(matrix.file, sep = ",", row.names=1) %>% as.matrix()
    #data <- data[,-ncol(data)]
    #data <- data %>% column_to_rownames("sample/gene") %>% as.matrix()
    data[data > 0] <- 1
    
    return(list(data, is.case))
}

permuteMatrix <- function(matrix.data, perm) {
    print(perm)
    data <- matrix.data[[1]]
    is.case <- matrix.data[[2]]
    n.samples <- length(is.case)
    n.cases <- length(which(is.case))
    n.ctrls <- n.samples - n.cases
    
    set.seed(perm)
    K <- sample.int(n.samples, size = n.cases, replace = FALSE)
    qv.cases <- rowSums(data[,K])
    qv.ctrls <- rowSums(data[,-K])
    
    gsum <- data.frame(qv.cases, n.cases - qv.cases, qv.ctrls, n.ctrls - qv.ctrls, stringsAsFactors = FALSE)
    gsum <- gsum %>% rownames_to_column(var = "Gene Name")
    colnames(gsum) <- c("Gene Name", "Qualified Case", "Unqualified Case", "Qualified Ctrl", "Unqualified Ctrl")
    return(gsum)
}


useCMHp <- function(gsum, gene) {
#    print(gene)
    conf <- lapply(gsum, function (x) x %>% filter(`Gene Name` == gene) %>% 
                       select(`Qualified Case`, `Unqualified Case`, `Qualified Ctrl`, `Unqualified Ctrl`))
    confa <- simplify2array(lapply(conf, function(x) array(unlist(x), dim = c(2,2))))
#    print(confa)

    cmh <- mantelhaen.test(confa, exact = TRUE)
    tcmh <- as_tibble(c(tidy(cmh)))
    
    return(tcmh)
}

doPermutations <- function(matrix.data, perm) {
    gsum <- lapply(matrix.data, function(x) permuteMatrix(x, perm))
    names(gsum) <- names(matrix.data)
    
    gene.names <- lapply(gsum, function(x) x$`Gene Name`)
    gene.names.r <- Reduce(union, gene.names)
    
    gsum <- lapply(gsum, function(x) x %>% full_join(as.tibble(gene.names.r), by = c("Gene Name" = "value")) %>% 
                       mutate(`Qualified Case` = case_when(is.na(`Qualified Case`) ~ 0, TRUE ~ as.double(`Qualified Case`)),
                              `Unqualified Case` = case_when(is.na(`Unqualified Case`) ~ max(x$`Qualified Case` + x$`Unqualified Case`), TRUE ~ (`Unqualified Case`)),
                              `Qualified Ctrl` = case_when(is.na(`Qualified Ctrl`) ~ 0, TRUE ~ as.double(`Qualified Ctrl`)),
                              `Unqualified Ctrl` = case_when(is.na(`Unqualified Ctrl`) ~ max(x$`Qualified Ctrl` + x$`Unqualified Ctrl`), TRUE ~ (`Unqualified Ctrl`))
                       ))    
    
    cmhs <- lapply(gene.names.r, function(i) useCMHp(gsum, i))
    
    cmhs.tidy <- unnest(tibble(cmhs)) %>% 
        mutate(gene = gene.names.r) %>% arrange(p.value) %>%  
        select(gene, p.value)
    
    permut_path <- paste0("Permutations_", gsub("_", "",  dir), "/")
    permut_path <- gsub("_/", "/", permut_path)
     
    saveRDS(cmhs.tidy, paste0(permut_path, dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS"))
}

getCMHresP <- function(model, dir, pc_path, nclust, nperm = c(1), cores = 1) {
    
    logfile <- paste0("Log/messagesperm_", dir, model, "_", nperm[1], "_", Sys.Date(), ".log")
    con <- file(logfile, "w")
    sink(con, append = TRUE, type = "message")
    print(model)
    
    matrix.data <- lapply(nclust, function(x) readMatrix(model, dir, pc_path, x))
    names(matrix.data) <- as.character(nclust)
    
    mclapply(nperm, function(x) doPermutations(matrix.data, x), mc.cores = cores)
    sink(type="message")
    close(con)
}
