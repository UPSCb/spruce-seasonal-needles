#' ---
#' title: "Principal Component Analysis (of a single calendar year)"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```
#' Libraries
suppressPackageStartupMessages(library(plotly))

###' PCA for vst blind counts
#' load samples infos and normalized (variance stabilized) gene expression data stored in rda file
#load("vst_blind.rda")
load("vst_aware.rda")

#' Focus on one year of sampling - samples past 2012-04-24 are the
#' next generation of needles
sel <- ordered_dates %in% levels(ordered_dates)[1:(grep("2012-05",levels(ordered_dates))[1]-1)]
ordered_dates <- factor(as.character(ordered_dates)[sel])
vst_aware <- vst_aware[,sel]

#' Compute principal component analysis
#pc <- prcomp(t(vst_blind)) 
pc <- prcomp(t(vst_aware)) 
# 't' transposes a matrix which means switch rows with columns
# before in vst.kt: rows=genes and columns=samples and data inside  = counts                    
# here PCA is calculated row after row, and we want to calculate it for each sample (so 84 PCs)

percent <- round(summary(pc)$importance[2,]*100)
#take the proportion of variance *100 and round for every PC
percent

#' plot according to PC1 and PC2  
p <- data.frame(sample = rownames(pc$x), pc$x[, 1:2], date = ordered_dates) %>% 
  ggplot(aes(x = PC1, y = PC2, colour = date)) +
  theme(legend.position = "right") + 
  labs(x=paste0("First component - ",percent[1],"%"),
       y=paste0("Second component - ", percent[2],"%")) + 
  geom_point()

ggplotly(p)

ggsave("pca_plot2_PC1_PC2_aware_1.png")

#' plot according to PC2 and PC3

p2 <- data.frame(sample = rownames(pc$x), pc$x[, 2:3], 
                 date = ordered_dates) %>% 
  ggplot(aes(x = PC2, y = PC3, colour = date)) +
  theme(legend.position = "right") +
  labs(x=paste0("Second component - ",percent[2],"%"),
       y=paste0("Third component - ", percent[3],"%")) + 
  geom_point()

ggplotly(p2)

ggsave("pca_plot2_PC2_PC3_aware_1.png")


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
