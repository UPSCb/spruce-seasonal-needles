#' ---
#' title: "Heatmap by season"
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
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")
cpal <- rainbow(30)
hpal <- colorRampPalette(colors = c("blue","white","red"))(100)

#' (Meta)Data
load("vst_aware.rda")
load("seasons-for-vst-aware.rda")

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

#' # Heatmap
#' ## All seasons
#' **Ordered by date**
lab <- as.character(ordered_dates)
lab[duplicated(lab)] <- ""
heatmap.2(t(scale(t(as.matrix(vst_aware[sel,])))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          Colv = FALSE,labCol = lab,
          dendrogram="row",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45)

season.lab <- as.character(seasons$season[match(sapply(strsplit(as.character(ordered_dates),"-"),"[",2),
                                   seasons$month)])
season.lab[duplicated(season.lab)] <- ""

heatmap.2(t(scale(t(as.matrix(vst_aware[sel,])))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          Colv = FALSE,labCol = season.lab,
          dendrogram="row",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45)

#' **Clustered**
heatmap.2(t(scale(t(as.matrix(vst_aware[sel,])))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          las=2,col=hpal,labCol = lab,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45)

#' ## By season
#' **Individual season, ordered by dates**
dev.null <- sapply(as.character(unique(seasons$season)),function(s,seasons){
  dat <- t(scale(t(as.matrix(vst_aware[sel,seasons==s]))))
  dat <- dat[rowSums(is.na(dat)) != ncol(dat),]
  heatmap.2(dat,distfun = pearson.dist,
            hclustfun = function(X){hclust(X,method = "ward.D2")},
            labRow=NA,trace="none",dendrogram = "row",
            Colv = FALSE,labCol = lab[seasons==s],
            key=TRUE,main=s,las=2,col=hpal,
            cexCol = 0.8,srtCol = 45)
},seasons$season[match(sapply(strsplit(as.character(ordered_dates),"-"),"[",2),
                     seasons$month)])

#' ## All seasons, average expression by season
seasonal.mean <- sapply(split.data.frame(t(vst_aware[sel,]),
                                         seasons$season[match(sapply(strsplit(as.character(ordered_dates),"-"),"[",2),
                                                                           seasons$month)]),colMeans)

seasonal.mean <- seasonal.mean[,c("early.summer","late.summer","autumn","winter","spring")]

#' **Column clustering**
heatmap.2(t(scale(t(seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          las=2,col=hpal,
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45)

#' **Sorted by season**
heatmap.2(t(scale(t(seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          las=2,col=hpal,
          Colv = FALSE, dendrogram="row",
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45)

pdf(file="expression_profiles/seasonal-average-expression-all-genes.pdf",width=24,height = 16)
heatmap.2(t(scale(t(seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=NA,trace="none",
          las=2,col=hpal,
          Colv = FALSE, dendrogram="row",
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45)
dev.off()

#' ## Photosynthesis genes
goi <- scan("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv", what="character")
goi_names <- scan("~/Git/UPSCb/projects/spruce-needles/doc/GOI_names.txt", what="character")

#' One gene is no part of the quantified genes - probably because it is classified as a putative contaminant
goi.sel <- which(goi %in% rownames(vst_aware))

#' One gene is never expressed
sdat <- t(scale(t(as.matrix(vst_aware[goi[goi.sel],]))))
na.sel <- which(is.na(rowSds(sdat)))
goi_names[goi.sel][which(is.na(rowSds(sdat)))]

#' **Ordered by date**
heatmap.2(sdat[-na.sel,],
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          Colv = FALSE,labCol = lab,
          dendrogram="row",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)

heatmap.2(sdat[-na.sel,],
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          Colv = FALSE,labCol = season.lab,
          dendrogram="row",
          las=2,col=hpal,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)

#' **Clustered**
heatmap.2(sdat[-na.sel,],
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          las=2,col=hpal,labCol = lab,
          ColSideColors = cpal[as.integer(factor(ordered_dates))],
          key=TRUE,main="Yearly expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)

#' **Seasonal average, sorted by season**
goi.seasonal.mean <- sapply(split.data.frame(t(vst_aware[goi[goi.sel],][-na.sel,]),
                                         seasons$season[match(sapply(strsplit(as.character(ordered_dates),"-"),"[",2),
                                                              seasons$month)]),colMeans)

goi.seasonal.mean <- goi.seasonal.mean[,c("early.summer","late.summer","autumn","winter","spring")]

heatmap.2(t(scale(t(goi.seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          las=2,col=hpal,
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)

heatmap.2(t(scale(t(goi.seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          las=2,col=hpal,
          Colv = FALSE, dendrogram="row",
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)

pdf(file="expression_profiles/seasonal-average-expression-photosynthetic-genes.pdf",width=24,height = 16)
heatmap.2(t(scale(t(goi.seasonal.mean))),
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method = "ward.D2")},
          labRow=goi_names[goi.sel][-na.sel],trace="none",
          las=2,col=hpal,
          Colv = FALSE, dendrogram="row",
          key=TRUE,main="Seasonal mean expression",
          cexCol = 0.8,srtCol = 45,cexRow = 0.6)
dev.off()

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
