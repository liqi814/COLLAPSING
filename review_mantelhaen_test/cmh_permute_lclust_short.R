# adapt to your path
setwd("/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/")
source("scripts/cmh_permute_lclust_short_functions.R")

dir <- "All"

pc_path <- "/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/flashpcaAll"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

permut_path <- "LatinoClustersAll/PermutationsAll"

if(!dir.exists(permut_path)) {
  dir.create(permut_path)
} else {
  print((paste0("Directory ", permut_path, "already exists")))
}

if(!dir.exists(paste0(permut_path, "/../Log"))) {
  dir.create(paste0(permut_path, "/../Log"))
} else {
  print((paste0("Directory Log already exists")))
}

LatinoIGMnclust <- c(3, 6) 
IGMnclust <-  c(0:10)

dbGAPnclust <-  c(0:6)
LatinodbGAPnclust <- c(3) 

IGM_resolution <- "0_25"
dbGAP_resolution <- "0_3"

#args = commandArgs(trailingOnly=TRUE)
#model <- args[1]
#model <- "dominantPTV"
#model <- "dominantSynonymous"
#model <- "dominantUltraRareP"
#model <- "dominantRareNB"
model <- "dominantRareEnsemble"
#model <- "dominantUltraRareEnsemble"

# you can start with 1:100 and later run 101:1000 if you want since the number is used as seed
nperm <- 1:250
cores <- 5

# run for each model separately
getCMHresP(model, dir, pc_path, permut_path, LatinoIGMnclust, LatinodbGAPnclust, IGM_resolution, dbGAP_resolution, nperm, cores)

# to check if there are permutations missing
files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(model, "_CMH_permut_perm_", "*.RDS")))
success <- as.numeric(gsub(paste0(model, "_CMH_permut_perm_|.RDS"), "", files))
nperm <- nperm[!nperm%in%success]
