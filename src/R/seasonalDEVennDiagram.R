#' ---
#' title: "Spruce needles seasonal differential expression analysis - figures"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/")

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/")
#' ```

#' Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(vsn))

#' Helpers
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/featureSelection.R"))
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/gopher.R"))
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' Palette
hpal <- colorRampPalette(colors = c("blue","red"))(100)
pal <- c("summer"="#48AF4D",
           "autumn"="#EC7063",
        "winter"="#D6EAF8",
          "spring"="#D7EB73")

#' # Differential Expression
#' ## Data
#' * Metadata
load(file="seasons-for-vst-aware.rda")
levels(seasons$season)[levels(seasons$season)=="late.summer"] <- "summer"

#' * Gene of Interest (GOI)
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv",
                header = FALSE,stringsAsFactors = FALSE)[,1]

#' * Raw data
load("kt_tech_merged_cropped_ordered.rda")

#' * Normalized data
load("vst_aware.rda")
#' Focus on one year of sampling - samples past 2012-04-24 are the
#' next generation of needles
sel <- ordered_dates %in% levels(ordered_dates)[1:(grep("2012-05",levels(ordered_dates))[1]-1)]
ord_dates <- factor(as.character(ordered_dates)[sel])
vst_aware <- vst_aware[,sel]

#' Reorder the samples
kt_tech_merged_cropped_ordered <- kt_tech_merged_cropped_ordered[,match(colnames(vst_aware),
                                                                        colnames(kt_tech_merged_cropped_ordered))]
stopifnot(colnames(vst_aware)==colnames(kt_tech_merged_cropped_ordered))

#' Remove the "early summer"
s.sel <- seasons$season[match(sapply(strsplit(as.character(ord_dates),"-"),"[",2),
                     seasons$month)] != "early.summer"

#' ## DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = kt_tech_merged_cropped_ordered[,s.sel],
                              colData = data.frame(season=seasons$season[match(sapply(strsplit(as.character(ord_dates),"-"),"[",2),
                                                                               seasons$month)][s.sel]),
                              design = ~season)

#' ## extend the colData
for (s in levels(colData(dds)$season)){
  colData(dds)[,s] <- factor(ifelse(colData(dds)$season==s,"yes","no"))
}

#' ## Differential Expression
dds <- DESeq(dds)

#' Cutoffs as devised from Schurch et al, RNA, 2016
#' We compare 6 (autumn) vs. 9 (winter) and 9 (winter) vs. 5 (spring), respectively
alpha=0.01
lfc=0.5

#' ### Every season vs. all other
resList <- mclapply(levels(colData(dds)$season),function(s){
  
  # change the design
  eval(parse(text=paste0("design(dds) <- ~",s)))
  
  # DE
  dds <- DESeq(dds)
  
  # results
  res <- results(dds)
  
  # MA plot
  DESeq2::plotMA(res)

  # Given what we know from the seasonal effect (overall less expressed genes in winter 
  # but which are on average more expressed due to the library size estimation), we should
  # consider larger fold changes as relevant. The seasonal effect analysis suggest a cutoff 
  # of at least 2 fold changes. The normalisation effect is actually visible in the volcano-plot
  # as the orange and yellow dots (at x=~0 and y=~0), that show a consistent shift from 0 towards
  # either season.
  volcanoPlot(res,lfc = 2)
  
  IDsUp <- rownames(res)[!is.na(res$padj) & res$padj <= alpha & res$log2FoldChange >= lfc]
  IDsDn <- rownames(res)[!is.na(res$padj) & res$padj <= alpha & res$log2FoldChange <= -lfc]
  IDs <- c(IDsUp,IDsDn)
  write(IDs,file=paste0("analysis/DESeq2/",s,"-vs-others_DE-gene-IDs_lfc-dot-five_FDR-one-percent.txt"))
  
  sel <- !is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= lfc
  write.table(data.frame(ID=rownames(res)[sel],log2FC=res[sel,"log2FoldChange"],FDR=res[sel,"padj"]), quote=FALSE,
              row.names = FALSE,
              sep = "\t",file=paste0("analysis/DESeq2/",s,"-vs-others_DE_lfc-dot-five_FDR-one-percent.txt"))
  
  sel <- !is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= lfc
  sel2 <- rownames(res)[sel] %in% goi
  write.table(data.frame(ID=rownames(res)[sel][sel2],log2FC=res[sel,"log2FoldChange"][sel2],FDR=res[sel,"padj"][sel2]), 
              quote=FALSE,
              row.names = FALSE,
              sep = "\t",file=paste0("analysis/DESeq2/",s,"-vs-others_DE_lfc-dot-five_FDR-one-percent.txt"))
  
  lfc=2
  IDsUp2 <- rownames(res)[!is.na(res$padj) & res$padj <= alpha & res$log2FoldChange >= lfc]
  IDsDn2 <- rownames(res)[!is.na(res$padj) & res$padj <= alpha & res$log2FoldChange <= -lfc]
  IDs2 <- c(IDsUp2,IDsDn2)
  write(IDs2,file=paste0("analysis/DESeq2/",s,"-vs-others_DE-gene-IDs_lfc-two_FDR-one-percent.txt"))
  sel <- !is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) >= lfc
  write.table(data.frame(ID=rownames(res)[sel],log2FC=res[sel,"log2FoldChange"],FDR=res[sel,"padj"]), quote=FALSE,
              row.names = FALSE,
              sep = "\t",file=paste0("analysis/DESeq2/",s,"-vs-others_DE_lfc-two_FDR-one-percent.txt"))

  list(IDsUp,IDsDn,IDs,IDsUp2,IDsDn2,IDs2)
  
},mc.cores=4L)

names(resList) <- levels(colData(dds)$season)

#' ## Plots
pdf(file="analysis/DESeq2/VennDiagrams.pdf")

#' ### VennDiagram all DE genes
#' #### The lfc cutoff is 0.5
#' 
#' * Up-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",1),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Up-regulated genes (lfc >= 0.5)"))

#' * Down-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",2),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Down-regulated genes (lfc <= -0.5)"))

#' * all DE
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",3),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="All genes (abs(lfc) >= 0.5)"))

#' #### The lfc cutoff is 2
#'
#' * Up-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",4),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Up-regulated genes (lfc >= 2)"))

#' * Down-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",5),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Down-regulated genes (lfc <= -2)"))

#' * all DE
grid.newpage()
grid.draw(venn.diagram(lapply(resList,"[[",6),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="All genes (abs(lfc) >= 2)"))

#' ### VennDiagram GOI only
#' #### The lfc cutoff is 0.5
#' * Up-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",1),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Up-regulated genes (lfc >= 0.5)"))

#' * Down-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",2),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Down-regulated genes (lfc <= -0.5)"))

#' * all DE
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",3),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="All genes (abs(lfc) >= 0.5)"))

#' #### The lfc cutoff is 2
#' * Up-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",4),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Up-regulated genes (lfc >= 2)"))

#' * Down-regulated only
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",5),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="Down-regulated genes (lfc <= -2)"))

#' * all DE
grid.newpage()
grid.draw(venn.diagram(lapply(lapply(resList,"[[",6),function(l){l[l %in% goi]}),
                       filename = NULL,
                       fill=pal[names(resList)],
                       category.names = names(resList),
                       main="All genes (abs(lfc) >= 2)"))
dev.off()

#' ## Enrichment
cutoffs <- sapply(seq(0,10,.5),featureSelect,counts=vst_aware[,s.sel],
                  conditions = colData(dds)$season,nrep=4)
plot(seq(0,10,.5),colSums(cutoffs),type="l",xlab="vst cutoff",ylab="# genes")

pop <- sub("\\.[0-9]+$","",rownames(vst_aware)[cutoffs[,2]])
  
#' ### Common genes
#' They have to do with cellulose synthase and laccase
enr <- gopher(genes=sub("\\.[0-9]+$","",Reduce(intersect,lapply(resList,"[[",3))),
              background = pop,url = "pabies")

enr$go[,c("id","padj")]
enr$go$name
enr$kegg$name
enr$pfam$name

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
