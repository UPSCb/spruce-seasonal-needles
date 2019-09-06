#' ---
#' title: "Seasonal heatmap clusters"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```

#' Libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyverse))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")
cpal <- rainbow(30)
hpal <- colorRampPalette(colors = c("blue","white","red"))(100)

#' (Meta)Data
load("vst_aware.rda")
load("seasons-for-vst-aware.rda")
samples <- read.csv(file="~/Git/UPSCb/projects/spruce-needles/doc/rna_seq_spruce_needles_2011_2012_sample_info.csv",
                    stringsAsFactors = FALSE)
samples <- samples[match(colnames(vst_aware),samples$Sampling_ID),]

#' # Pre-process
#' ## Filter
#' Focus on one year of sampling - samples past 2012-04-24 are the
#' next generation of needles
sel <- ordered_dates %in% levels(ordered_dates)[1:(grep("2012-05",levels(ordered_dates))[1]-1)]
ordered_dates <- factor(as.character(ordered_dates)[sel])
vst_aware <- vst_aware[,sel]

#' ## Expression cutoff
#' Select genes that are expressed above a treshold of counts 'exp' and in at least 'n' replicates
plot(sapply(X = 1:10, function(i){sum(featureSelect(vst_aware,ordered_dates,i))}),type="l",
     main="number of genes at different threshold of expression",xlab="threshold",ylab="number of genes")

#' Selecting genes whose counts are over or equal to a vst value of 1.5 in at least 2 replicates of one time point
sel <- featureSelect(vst_aware,ordered_dates,1.5)

#' # Plots
#' ## All seasons, average expression by season
seasonal.mean <- sapply(split.data.frame(t(vst_aware[sel,]),
                                         seasons$season[match(sapply(strsplit(as.character(ordered_dates),"-"),"[",2),
                                                              seasons$month)]),colMeans)

seasonal.mean <- seasonal.mean[,c("early.summer","late.summer","autumn","winter","spring")]

hc <- hclust(pearson.dist(t(scale(t(seasonal.mean)))),"ward.D2")
test <- cutree(hc,12)
barplot(table(test))
library(LSD)
pdf("twelve-clusters.pdf",width=12,height=8)
par(mfrow=c(2,1))
dev.null <- sapply(1:12,function(i){
  g <- names(test[test==i])

  test2 <- split.data.frame(scale(t(vst_aware[sel,][g,])),samples$Sampling_Date)
  stopifnot(all(names(test2) == levels(ordered_dates)))
  
  
  
  DF <- data.frame(time=as.Date(levels(ordered_dates)),
                   t(sapply(test2,quantile)))
  
  p <- ggplot(DF, aes(time, group = 1)) +
    annotate("rect", xmin=as.Date("2011-10-01"),
             xmax=as.Date("2012-03-01"),
             ymin=-Inf,ymax=+Inf,
             fill="#ADD8E6",alpha=0.3) +
    geom_line(aes(y=X50.), color="blue") +
    geom_ribbon(aes(ymin=X0., ymax=X100.), alpha=0.3) + 
    geom_ribbon(aes(ymin=X25., ymax=X75.), alpha=0.6) + 
    scale_x_date(date_breaks = "month") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(y = "z-score") + ggtitle(paste("Expression profile of cluster", i),
                                     subtitle = sprintf("Average vst expression: %s - of %s genes",
                                                      round(mean(vst_aware[sel,][g,]),digits=2),
                                                      length(g)))
  
  plot(p)
  
  #suppressMessages(ggsave(filename = paste0("expression_profiles/",goi_names[i,1],"_expression_profile.jpeg")))
  
  })
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
