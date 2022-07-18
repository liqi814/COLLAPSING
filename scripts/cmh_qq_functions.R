library(tidyverse)
# library(dplyr)
# library(tidyr)
# library(tibble)
# library(readr)
# library(ggplot2)
# library(purrr)

library(broom)
library(readxl)
library(data.table)
library(ggrepel)

readPermTable <- function(dir, model, perm) {
    print(perm)

    perm_file <- paste0(permut_path, dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")

     if(!file.exists(perm_file)) {
        if(!str_detect(model, "recessive|dominant")) {
            model <- paste0("dominant", model)
            perm_file <- paste0(permut_path, dir, model, "_CMH_permut_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")
        }
        if(!file.exists(perm_file)) {
            perm_file <- paste0(permut_path, dir, model, "_CMH_perm_res_", resolution, "_min_sample_", min_sample, "_perm_", perm, ".RDS")
            if(!file.exists(perm_file)) {
            perm_file <- paste0(permut_path, dir, model, "_CMH_perm_res_", resolution, "_genderstrat_min_sample_", min_sample, "_perm_", perm, ".RDS")            
            }
        }
    }

    if(file.exists(perm_file)) {
        cmhp <- readRDS(perm_file)
        return(cmhp$p.value)
    } else {
        print(perm_file)
        return(NULL)
    }

}

plotQQ <- function(dir, model, cmh, nperm) {
    permRes <- lapply(nperm, function(x) readPermTable(dir, model, x))
    permRes.matrix <- do.call(cbind, permRes)
    permRes.matrix <- permRes.matrix[1:nrow(cmh),]
    pVal <- rowMeans(permRes.matrix)
    quant <- apply(permRes.matrix, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
    
    obs <- cmh$p.value
    exp <- pVal
    
    obs <- obs[which(!is.na(obs))]
    exp <- exp[which(!is.na(exp))]
    obs[obs > 1] <- 1
    exp[exp > 1] <- 1
    exp <- exp[1:length(obs)]
    
    gws <- 0.05 / length(obs)
    reg_obs <- obs[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_exp <- exp[(obs > gws) & (exp > gws) & (exp < 1) & (obs < 1)]
    reg_obs <- qchisq(reg_obs, 1, lower.tail = FALSE)
    reg_exp <- qchisq(reg_exp, 1, lower.tail = FALSE)
    
    obs <- -log10(obs)
    exp <- -log10(exp)
    
    uPerc <- -log10(quant[1,])[1:length(obs)]
    lPerc <- -log10(quant[2,])[1:length(obs)]
    
    d <- data.frame(gsub("'", "", cmh$gene)[!is.na(cmh$p.value)], obs, exp, uPerc, lPerc)
    names(d) <- c("gene", "obs", "exp", "upper", "lower")
    
    lin.reg <- lm(reg_obs ~ 0 + reg_exp[1:length(reg_obs)])
    lambda <- lin.reg$coefficients
   
    save(list=ls(), file=paste0(qq_path, model, "_qq_plot_res_", resolution, "_min_sample_", min_sample, ".RData"))
 
    ggplot(d, aes(exp, obs)) +
        geom_point(colour="red") +
        geom_text_repel(data = subset(d, obs > 6), aes(label = gene)) +
        geom_line(aes(exp, upper), colour="orange") +
        geom_line(aes(exp, lower), colour="green") +
        geom_abline(slope = 1, colour="blue") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(paste0("QQ Plot: Observed vs. expected p-values. Lambda = ", round(lambda, digits = 4))) +
        xlab("Expected -log10(p)") +
        ylab("Observed -log10(p)")
    
    ggsave(filename = paste0(qq_path, dir, model, "_CMH_qq_exact_res_", resolution, "_min_sample_", min_sample, ".pdf"), 
           device = "pdf", width = 8, height = 8)
    
}        

plotQQAll <- function(dir, model, cmhs, nperm) {
    print(model)
	# tab names in xslx files can only be 30 characters
    if (model == "dominantFlexiblePolyphenDamaging") {
        cmh <- cmhs[["dominantFlexiblePolyphenDamagin"]]
    } else {
        cmh <- cmhs[[model]]
    }
    plotQQ(dir, model, cmh, nperm)
}
