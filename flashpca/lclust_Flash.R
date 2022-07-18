#library(tidyverse)
library(tidyr)
library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(umap)
library(RColorBrewer)


#######adapt paths to your structure
pc_path <- "/nfs/projects/dbGap/COPD_sample_remove/flashpca"

evf <- list.files(pc_path, "*flashpca_eigenvectors$", full.names = TRUE)[1]
pcs <- fread(evf)
pcs <- pcs %>% rename_at(vars(U1:U10), ~ paste0("PC", 1:10))

sf <- "/nfs/projects/dbGap/COPD_sample_remove/KING/atav_prune/dbGap_noIGM_noCOPD_pruned_sample.txt"
samples <- fread(sf, header = FALSE)

cases <- samples %>% filter(V6 == 2)
controls <- samples %>% filter(V6 == 1)

pcs <- pcs %>% left_join(samples, by = c("IID" = "V2")) %>% 
    mutate(Phenotype = case_when(IID %in% cases$V2 ~ "case",
                                 IID %in% controls$V2 ~ "control",
                                 TRUE ~ "other"),
           Gender = case_when(V5 == 1 ~ "male",
                              V5 == 2 ~ "female",
                              TRUE ~ "other")) %>% 
    select(-c(V1:V8))

#######defaults##########
args<-commandArgs(TRUE)
resolution <- as.numeric(args[1])  # 0.1-0.4
print(resolution)
typeof(resolution)
resolution2 <- gsub("\\.", "_", as.character(resolution))

dir <- "All"

if (dir == "") {
    dir2 = ""
} else {
    dir2 <- paste0(dir, "_", resolution2, "_")
}


plot_path <- paste0("Plots_", dir, "/")
plot_path <- gsub("_/", "/", plot_path)

if(!dir.exists(plot_path)) {
    dir.create(plot_path)
} else {
    print((paste0("Directory ", plot_path, "already exists")))
}


#### Louvain Clustering ####
total <- pcs %>% select(PC1:PC6) %>% as.matrix()
rownames(total) <- pcs$IID

total_nn <- FindNeighbors(total)
total_snn <- Matrix(total_nn[[2]], sparse=TRUE)

total_clust <- FindClusters(total_snn, resolution = resolution)

saveRDS(total_clust, paste0(plot_path, dir2, "_lclust_res_", resolution2, ".RDS"))

# total_clust <- readRDS(paste0(plot_path, dir2, "lclust_res_", resolution2, ".RDS"))

pcs <- pcs %>% mutate(cluster = total_clust[,1])

#### visualize with UMAP ####
custom.settings <- umap.defaults
custom.settings$min_dist <- 0.1
custom.settings$random_state <- 188 # try another seed if the plot is not nice
custom.settings$verbose <- TRUE
custom.settings$metric <- "pearson"

umap <- umap(total, config = custom.settings)
saveRDS(umap, paste0(plot_path, dir2, "umap.RDS"))

y <- data.frame(umap$layout)
y <- y %>% mutate(clusters = total_clust[,1])

ncolors <- length(unique(pcs$cluster))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(ncolors)

ggplot(y, aes(x = X1, y = X2, color = clusters)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() + 
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
ggsave(filename = paste0(plot_path, dir2, "UMAP_", resolution2, ".png"), 
       device = "png", width = 6, height = 5)

ggplot(y, aes(x = X1, y = X2, color = clusters)) +
    geom_point(size = 0.5, alpha = 0.6) +
    theme(legend.title=element_blank()) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    theme_minimal() +
    scale_color_manual(values = mycolors) +
    guides(colour = guide_legend(override.aes = list(size = 3)))+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA)
        )
ggsave(filename = paste0(plot_path, dir2, "UMAP_", resolution2, "_transparent_bg.pdf"),
       device = "pdf", width = 6, height = 5,bg = "transparent")

#### table of cluster sizes ####
cl_ca_co <- pcs %>% group_by(cluster, Phenotype) %>% tally() %>% spread(Phenotype, n)
write.table(cl_ca_co, paste0(plot_path, dir2, "lclustering_res_", resolution2, "_cluster_sizes.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


#### write out cluster sample files ####
for (cl in levels(pcs$cluster)) {
    ind <- pcs %>% filter(cluster == cl)
    samples_cl <- samples %>% filter(V2 %in% ind$IID)
    write.table(samples_cl, sep = "\t", quote=FALSE, row.names = FALSE, col.names = FALSE, 
                file=paste0(pc_path, "/flashPCA_lclustering_res_", resolution2, "_cluster_", cl, "_sample.txt"))
    
}
