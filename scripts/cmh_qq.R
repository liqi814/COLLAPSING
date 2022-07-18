# adapt path
setwd("/nfs/projects/dbGap/COPD_sample_remove/collapsing_analyses/")
source("scripts/cmh_qq_functions.R")


dir <- "All"

if (dir == "") {
    dir = ""
} else {
    dir <- paste0(dir, "_")
}

plot_path <- paste0("/nfs/projects/dbGap/COPD_sample_remove/flashpca/Plots_", gsub("_", "",  dir), "/")
plot_path <- gsub("_/", "/", plot_path)
results_path <- paste0("CMH_", gsub("_", "",  dir), "/")
results_path <- gsub("_/", "/", results_path)
qq_path <- paste0("QQ_plot_", gsub("_", "",  dir), "/")
qq_path <- gsub("_/", "/", qq_path)
permut_path <- paste0("Permutations_", gsub("_", "",  dir), "/")
permut_path <- gsub("_/", "/", permut_path)

if(!dir.exists(qq_path)) {
        dir.create(qq_path)
    } else {
        print((paste0("Directory ", qq_path, " already exists")))
    }


##setting
resolution <- "0_3"
min_sample <- 5
args = commandArgs(trailingOnly=TRUE)
model <- args[1]
#model <- "RareEnsemble"


# adapt to your filename
cl_sizes <- fread(paste0(plot_path, dir, resolution, "_lclustering_res_", resolution, "_cluster_sizes.txt"))
nclust <- (cl_sizes %>% filter(case >= min_sample & control >= min_sample))$cluster

cmh <- read_xlsx(paste0(results_path, dir, model,"_CMH_exact_summary_lclust_res_", resolution,".xlsx"))

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
