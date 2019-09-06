#' ---
#' title: "Soile's spruce seasonal wood series manuscript heatmaps"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/htuominen/27_SpruceTimeSeries")
#' ```

#' Libraries
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(RColorBrewer))

#' palette
# pal <- rev(brewer.pal(9,"RdBu"))
pal <- colorRampPalette(colors = c("blue","white","red"))(100)
cpal <- c(phloem=rgb(27,158,119,maxColorValue=255),xylem=rgb(255,0,255,maxColorValue=255))

#' Read the data
#load("analysis/HTSeq/20161216/library-size-normalized_variance-stabilized_data.rda")
#load("analysis/kallisto/raw-unormalised-gene-expression_data.csv")
load("analysis/kallisto/combined-library-size-normalized_variance-stabilized_gene-expression_data.csv")

#' Read the doc
#samples <- read.csv("~/Git/UPSCb/projects/spruce-wood-time-series/doc/samplesmod.csv",row.names=1)
samples <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/rna_seq_spruce_needles_2011_2012_sample_info.csv")

#' Order the Date factor by date and not alphabetically
# dates <- rev(unique(as.character(samples$Sampling_Date)))
dates <- unique(as.character(samples$Sampling_Date))
#dates <- dates[c(3:23,1,2)]
# for(r in dates){
#   samples$Sampling_Date <- relevel(samples$Sampling_Date,r)
# }

#' Sort by Date and then by Tissue
# ord <- order(samples$Date,samples$Tissue)
# samples <- samples[ord,]
# vst <- vst[,ord]
# stopifnot(all(rownames(samples) == colnames(vst)))

#' Create a label
lbls <- as.character(samples$Sampling_Date)
lbls[duplicated(lbls)] <- NA
#lbls <- sub("-201[1,2]$","",lbls)

#' # Process
#' Calculate the standard deviation per gene per sample replicates
sds <- do.call("cbind",
               mclapply(
                 split.data.frame(
                   t(vst.kt),
                   #paste0(samples$Tissue,samples$Sampling)
                   ),
                 colSds,
                 mc.cores=4))
sds[is.na(sds)] <- 0

#' Select those genes which expression is constant across all replicates
sel <- rowSums(sds <= 1) == 44 & (phloem | xylem)

#' ## Saturated expression
vst.sat <- vst[sel,]
vst.sat[vst.sat > 10 ] <- 10
vst.sat[vst.sat < 3 ] <- 3

#' # Plot
#' ## By date (main)
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          Colv = FALSE, dendrogram = "row",
          las=2,col=pal,
          ColSideColors = cpal[as.integer(samples$Tissue)],
          labCol = lbls,
          key=FALSE,main="Saturated expression (3 to 10, all samples)",
          cexCol = 0.8,srtCol = 45)

#' Heatmap of all replicates across one season. Samples are ordered by sampling date, the 
#' expression values (variance stabilised, library sized normalised, i.e. approx. log2)
#' were saturated to enhance the legibility of the heatmap. Expression values lower than 3 
#' and higher than 10 were raised or reduced to 3 or 10, respectively. The ribbon at the top
#' of the heatmap indicate the tissue type, magenta for xylem and green for phloem. 
#' Almost every time point consists of 3 replicates of both tissues, but only the first is
#' labelled by its sampling date.
#' 
#' The heatmap clearly shows that the sequencing data is of high quality, all replicates 
#' showing very similar patterns. It also shows that the data is adequately picking up
#' both seasonal and tissue trends. Most strikingly, is the shift of expression of 
#' many genes in the winter time, with surprisingly many genes getting a significantly higher
#' relative expression. It is also interesting to observe relative increase/decrease pattern
#' of expression of some gene clusters. Secondly, tissue specific expression can be observed
#' in the checkered pattern of the heatmap. There it is interesting to see that this patterns 
#' disappear during the winter. After winter, in early spring (lower part of the heatmap, 
#' sampling date in March and April), it would seem that the phloem gets re-activated earlier
#' than the xylem.
#' 
#' Finally, it is intersting to see in the supplementary heatmap figure SXX, where the samples were 
#' also clustered by column, that 4 high level clusters are evident. This correspond to not 4 but 3
#' "seasons". The growing (summer), the dormancy (winter), and the transition (spring and autumn) seasons. 
#' It is during the growing season, that the different tissues show the highest divergence in expression.
#' On the other hand, the winter, past october shows a very consistant expression. The term "dormancy" 
#' could even be challenged at the gene expression level, since as many genes are expressed as in the other
#' time of the year, with one large cluster being very specific to winter. Finally, it is interesting
#' to see the high similarity in patterns of expression during the transition period although September
#' (Autumn) separates from spring (March-May), it shows striking similarties. Nonetheless, effectively, all 4 seasons 
#' can be observed. 
#' 
#' Out of all samples, 2 seem to be slight outliers, both Phloem, one from June 27th and one from August 1st.


pdf(file="analysis/heatmap-by-date.pdf",paper = "a4r")
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          Colv = FALSE, dendrogram = "row",
          las=2,col=pal,
          ColSideColors = cpal[as.integer(samples$Tissue)],
          labCol = lbls,
          key=FALSE,main="Saturated expression (3 to 10, all samples)",
          cexCol = 0.8,srtCol = 45)
dev.off()

#' ## Clustered by Column (supp)
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          las=2,col=pal,
          ColSideColors = cpal[as.integer(samples$Tissue)],
          labCol = sub("-201[1,2]$","",as.character(samples$Date)),
          key=FALSE,main="Saturated expression (3 to 10, all samples)",
          cexCol = 0.6)

pdf(file="analysis/heatmap-supplementary.pdf",paper = "a4r")
heatmap.2(as.matrix(vst.sat),
          labRow=NA,trace="none",
          las=2,col=pal,
          ColSideColors = cpal[as.integer(samples$Tissue)],
          labCol = sub("-201[1,2]$","",as.character(samples$Date)),
          key=FALSE,main="Saturated expression (3 to 10, all samples)",
          cexCol = 0.6)
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
