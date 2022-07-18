# adapt to your path
setwd("/nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/")
source("scripts/cmh_permute_lclust_cohort_functions.R")

dir <- "All"

pc_path <- "/nfs/projects/dbGap/mantelhaen_test/review_permutations_EuropeanVSnon-Euro/flashpcaAll"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

permut_path <- "EuropeanClustersAll/PermutationsdbGAP"

if(!dir.exists(permut_path)) {
  dir.create(permut_path)
} else {
  print((paste0("Directory ", permut_path, " already exists")))
}

if(!dir.exists(paste0(permut_path, "/../Log"))) {
  dir.create(paste0(permut_path, "/../Log"))
} else {
  print((paste0("Directory Log already exists")))
}

EuroIGMnclust <- c(0, 1, 2, 10) # European clusters
IGMnclust <-  c(0:10)
NonEuroIGMnclust <- IGMnclust[ ! IGMnclust %in% EuroIGMnclust]

dbGAPnclust <-  c(0:6)
EurodbGAPnclust <- c(0, 1, 4, 6) # European clusters
NonEurodbGAPnclust <- dbGAPnclust [ ! dbGAPnclust %in% EurodbGAPnclust ]

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
getCMHresP(model, dir, pc_path, NonEuroIGMnclust, IGM_resolution, nperm, cores)

# to check if there are permutations missing
files <- list.files(paste0(permut_path), pattern = glob2rx(paste0(model, "_CMH_permut_perm_", "*.RDS")))
success <- as.numeric(gsub(paste0(model, "_CMH_permut_perm_|.RDS"), "", files))
nperm <- nperm[!nperm%in%success]
