#' ---
#' title: "Spruce needles seasonal differential expression analysis"
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
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(wordcloud))

#' Helpers
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/featureSelection.R"))
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/gopher.R"))
suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' Palette
hpal <- colorRampPalette(colors = c("blue","red"))(100)

#' # Pre-process
#' ## Seasons
#' Read (Meta)Data
load("vst_aware.rda")

#' ### Filter
#' Focus on one year of sampling - samples past 2012-04-24 are the
#' next generation of needles
sel <- ordered_dates %in% levels(ordered_dates)[1:(grep("2012-05",levels(ordered_dates))[1]-1)]
ord_dates <- factor(as.character(ordered_dates)[sel])
vst_aware <- vst_aware[,sel]

#' ### Hierarchical clustering
#' We define 5 seasons: 
#' 
#' * spring: March - April
#' 
#' * early-summer: May - June
#' 
#' * summer: July-August
#' 
#' * autumn: September
#' 
#' * winter: October-February
#' 
hc <- hclust(dist(t(vst_aware)))
plot(hc, main = "Hierarchical clustering",labels = ord_dates,cex=0.5)

seasons <- data.frame(month=sprintf("%02d",1:12),
                      season=c(rep("winter",2),
                               rep("spring",2),
                               rep("early.summer",2),
                               rep("late.summer",2),
                               "autumn",
                               rep("winter",3)))

save(seasons,file="seasons-for-vst-aware.rda")

#' ## Raw data
#' ### Raw counts
load("kt_tech_merged_cropped_ordered.rda")
kt_tech_merged_cropped_ordered <- kt_tech_merged_cropped_ordered[,match(colnames(vst_aware),colnames(kt_tech_merged_cropped_ordered))]
stopifnot(colnames(vst_aware)==colnames(kt_tech_merged_cropped_ordered))

#' ### DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = kt_tech_merged_cropped_ordered,
                              colData = data.frame(season=seasons$season[match(sapply(strsplit(as.character(ord_dates),"-"),"[",2),
                                                                               seasons$month)]),
                              design = ~season)

#' ## Differential Expression
dds <- DESeq(dds)

#' Cutoffs as devised from Schurch et al, RNA, 2016
#' We compare 6 (autumn) vs. 9 (winter) and 9 (winter) vs. 5 (spring), respectively
alpha=0.01
lfc=0.5

#' ### Autumn vs. Winter
res <- results(dds,c("season","autumn","winter"))

#' *MA plot*
DESeq2::plotMA(res)

#' *Volcanoplot*
#' 
#' Given what we know from the seasonal effect (overall less expressed genes in winter 
#' but which are on average more expressed due to the library size estimation), we should
#' consider larger fold changes as relevant. The seasonal effect analysis suggest a cutoff 
#' of at least 2 fold changes. The normalisation effect is actually visible in the volcano-plot
#' as the orange and yellow dots (at x=~0 and y=~0), that show a consistent shift from 0 towards
#' either season.
volcanoPlot(res)
lfc=2

#' *Tables*
dir.create("analysis/DESeq2",showWarnings = FALSE)
write.table(res,file="analysis/DESeq2/autumn-vs-winter_DE.txt",
            sep="\t",quote = FALSE)

avsw <- sub("\\.1$","",rownames(res)[!is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) > lfc])
write.table(avsw,file="analysis/DESeq2/autumn-vs-winter_DE-gene-IDs.txt",
            sep="\t",quote = FALSE,col.names = FALSE, row.names = FALSE)

#' ### Winter vs. Spring
res <- results(dds,c("season","winter","spring"))

#' *MA plot*
DESeq2::plotMA(res)

#' *Volcanoplot*
#' 
#' Given what we know from the seasonal effect (overall less expressed genes in winter 
#' but which are on average more expressed due to the library size estimation), we should
#' consider larger fold changes as relevant. The seasonal effect analysis suggest a cutoff 
#' of at least 2 fold changes. The normalisation effect is actually visible in the volcano-plot
#' as the orange and yellow dots (at x=~0 and y=~0), that show a consistent shift from 0 towards
#' either season.
volcanoPlot(res)
alpha=0.01
lfc=2

#' *Tables*
write.table(res,file="analysis/DESeq2/winter-vs-spring_DE.txt",
            sep="\t",quote = FALSE)

wvss <- sub("\\.1$","",rownames(res)[!is.na(res$padj) & res$padj <= alpha & abs(res$log2FoldChange) > lfc])
write.table(wvss,file="analysis/DESeq2/winter-vs-spring_DE-gene-IDs.txt",
            sep="\t",quote = FALSE,col.names = FALSE, row.names = FALSE)


#' ## Gene Ontology
#' ### The population
#' We select genes that are expressed above a vst value of in at least 2 replicates of any season
vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
#' The VST is satisfying
meanSdPlot(vst[rowSums(vst)>0,])

#' Devise the population
pop <- sub("\\.1$","",rownames(vst)[featureSelect(vst,colData(dds)$season)])
stopifnot(all(avsw %in% pop))
stopifnot(all(wvss %in% pop))

#' ### Autumn vs. Winter
enrichment <- gopher(avsw,
                     task = list("go","kegg","pfam"),
                     background = pop,url="pabies")

#' #### GO
#' Filter for OBSOLETE
enrichment$go <- enrichment$go[!grepl("^OBSOLETE.",enrichment$go$def),]
stopifnot(all(0 <= as.numeric(enrichment$go$padj) & as.numeric(enrichment$go$padj) <= 1))

#' Save
write.table(enrichment$go,file="analysis/DESeq2/autumn-vs-winter_GO.txt",
      sep="\t",quote = FALSE, row.names=FALSE)

write.table(enrichment$go[,c("id","padj")],file="analysis/DESeq2/autumn-vs-winter_REVIGO-input.txt",
            sep="\t",quote = FALSE, row.names=FALSE, col.names = FALSE)

#' #### KEGG
alpha=.1
write.table(enrichment$kegg,file="analysis/DESeq2/autumn-vs-winter_KEGG.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

#' Find the pathway of the significantly enriched enzymes
pathways <- keggLink("pathway",enrichment$kegg$id[enrichment$kegg$padj<=alpha])

#' Pathways are duplicated, clean up
stopifnot(length(grep("map",pathways)) == length(grep("ec",pathways)))
pathways <- pathways[grep("map",pathways)]

#' Proportion of the pathways
barplot(sort(table(sub("path:","",pathways)),decreasing = TRUE),las=2)

#' Get the pathway info
pinfo <- lapply(split(unique(pathways), ceiling(seq_along(unique(pathways))/10)),keggGet)

pathway.df <- do.call(rbind,lapply(pinfo,function(x){data.frame(ENTRY=sapply(x,"[[","ENTRY"),
                                                                NAME=sapply(x,"[[","NAME"))}))
tab <- table(sub("path:","",pathways))

pathway.df$OCCURENCE <- tab[match(pathway.df$ENTRY,names(tab))]

write.table(pathway.df,file="analysis/DESeq2/autumn-vs-winter_KEGG-pathway.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

#' Create a wordcloud of the pathway names
wordcloud(pathway.df$NAME,pathway.df$OCCURENCE/sum(pathway.df$OCCURENCE),
          colors = hpal,scale = c(2,.5),rot.per = 0)

#' #### PFAM
alpha <- 0.01
write.table(enrichment$pfam,file="analysis/DESeq2/autumn-vs-winter_PFAM.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

wordcloud(enrichment$pfam[enrichment$pfam$padj <= alpha,"name"],
          (as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"nt"])/
            as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"mt"])) * 
            (as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"nt"]) /
            as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"n"])),
          colors = hpal,scale = c(2,.5),rot.per = 0)

#' ### Winter vs. Spring
enrichment <- gopher(wvss,
                     task = list("go","kegg","pfam"),
                     background = pop,url="pabies")

#' #### GO
#' Filter out OBSOLETE terms
enrichment$go <- enrichment$go[!grepl("^OBSOLETE.",enrichment$go$def),]
stopifnot(all(0 <= as.numeric(enrichment$go$padj) & as.numeric(enrichment$go$padj) <= 1))

#' Save the results
write.table(enrichment$go,file="analysis/DESeq2/winter-vs-spring_GO.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

write.table(enrichment$go[,c("id","padj")],file="analysis/DESeq2/winter-vs-spring_REVIGO-input.txt",
            sep="\t",quote = FALSE, row.names=FALSE, col.names = FALSE)

#' #### KEGG
alpha=.1
write.table(enrichment$kegg,file="analysis/DESeq2/winter-vs-spring_KEGG.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

#' Find the pathway of the significantly enriched enzymes
pathways <- keggLink("pathway",enrichment$kegg$id[enrichment$kegg$padj<=alpha])

#' Pathways are duplicated, clean up
stopifnot(length(grep("map",pathways)) == length(grep("ec",pathways)))
pathways <- pathways[grep("map",pathways)]

#' Proportion of the pathways
barplot(sort(table(sub("path:","",pathways)),decreasing = TRUE),las=2)

#' Get the pathway info
pinfo <- lapply(split(unique(pathways), ceiling(seq_along(unique(pathways))/10)),keggGet)

pathway.df <- do.call(rbind,lapply(pinfo,function(x){data.frame(ENTRY=sapply(x,"[[","ENTRY"),
                                                          NAME=sapply(x,"[[","NAME"))}))

tab <- table(sub("path:","",pathways))

pathway.df$OCCURENCE <- tab[match(pathway.df$ENTRY,names(tab))]

write.table(pathway.df,file="analysis/DESeq2/winter-vs-spring_KEGG-pathway.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

#' Create a wordcloud of the pathway names
wordcloud(pathway.df$NAME,pathway.df$OCCURENCE/sum(pathway.df$OCCURENCE),
          colors = hpal,rot.per = 0)

#' #### PFAM
alpha <- 0.01
write.table(enrichment$pfam,file="analysis/DESeq2/winter-vs-spring_PFAM.txt",
            sep="\t",quote = FALSE, row.names=FALSE)

wordcloud(enrichment$pfam[enrichment$pfam$padj <= alpha,"name"],
          (as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"nt"])/
             as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"mt"])) * 
            (as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"nt"]) /
               as.integer(enrichment$pfam[enrichment$pfam$padj <= alpha,"n"])),
          colors = hpal,scale = c(2,.5),rot.per = 0)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```



