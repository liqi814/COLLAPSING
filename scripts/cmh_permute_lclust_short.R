# adapt to your path
setwd("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/")
source("scripts/cmh_permute_lclust_short_functions.R")

dir <- "All"

pc_path <- "/nfs/projects/dbGap/COPD_sample_remove/flashpca"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}


plot_path <- paste0("/nfs/projects/dbGap/COPD_sample_remove/flashpca/Plots_", gsub("_", "",  dir), "/")
plot_path <- gsub("_/", "/", plot_path)
#/nfs/projects/dbGap/COPD_sample_remove/flashpca/Plots_All/All_0_3_lclustering_res_0_3_cluster_sizes.txt
results_path <- paste0("CMH_", gsub("_", "",  dir), "/")
results_path <- gsub("_/", "/", results_path)
permut_path <- paste0("Permutations_", gsub("_", "",  dir), "/")
permut_path <- gsub("_/", "/", permut_path)


if(!dir.exists(permut_path)) {
  dir.create(permut_path)
} else {
  print((paste0("Directory ", permut_path, "already exists")))
}

if(!dir.exists("Log")) {
  dir.create("Log")
} else {
  print((paste0("Directory Log already exists")))
}

resolution <- "0_3"
min_sample <- 5

# adapt to your filename
cl_sizes <- fread(paste0(plot_path, dir, resolution, "_lclustering_res_", resolution, "_cluster_sizes.txt"))
nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster

args = commandArgs(trailingOnly=TRUE)
model <- args[1]
#model <- "dominantPTV"
#model <- "dominantSynonymous"
#model <- "dominantUltraRareP"
#model <- "dominantRareNB"
#model <- "dominantRareEnsemble"
#model <- "dominantUltraRareEnsemble"

# you can start with 1:100 and later run 101:1000 if you want since the number is used as seed
nperm <- 1:250
cores <- 5

# run for each model separately
getCMHresP(model, dir, pc_path, nclust, nperm, cores)

# to check if there are permutations missing
files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", "*.RDS")))
success <- as.numeric(gsub(paste0(dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_|.RDS"), "", files))
nperm <- nperm[!nperm%in%success]
