# adapt path
setwd("/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/LatinoClustersAll/")
source("../scripts/cmh_qq_functions.R")


dir <- "All"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

results_path <- "/nfs/projects/dbGap/mantelhaen_test/review_permutations_Latino/CMH_All/"
qq_path <- "QQplotAll/"
permut_path <- "PermutationsAll/"

if(!dir.exists(qq_path)) {
        dir.create(qq_path)
    } else {
        print((paste0("Directory ", qq_path, " already exists")))
    }



#args = commandArgs(trailingOnly=TRUE)
#model <- args[1]
#model <- "dominantPTV"
#model <- "dominantSynonymous"
#model <- "dominantUltraRareP"
#model <- "dominantRareNB"
model <- "dominantRareEnsemble"
#model <- "dominantUltraRareEnsemble"


# adapt to your filename

cmh <- read_xlsx(paste0(results_path, "All_Latino_", model, "_CMH_exact_summary_lclust.xlsx"))

#cmhs <- xlsx %>%
#    excel_sheets() %>%
#    purrr::set_names() %>%
#    map(read_excel, path = xlsx)

#models <- names(cmhs)
# tab names in xslx files can only be 30 characters
# models[grepl("dominantFlexiblePolyphen", models)] <- "dominantFlexiblePolyphenDamaging"

nperm <- 1:250

#lapply(models, function(x) plotQQAll(dir, x, cmhs, nperm))
plotQQ(dir, model, cmh, nperm)


# model <- "dominantUltraRareEnsemble"
# cmh <- readRDS(paste0(dir, model, "_res_", resolution, "_exact_CMH.RDS"))
# plotQQ(dir, model, cmh, nperm)
