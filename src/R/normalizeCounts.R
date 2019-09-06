#' ---
#' title: "Spruce Needles Kallisto Biological QA"
#' author: "Nicolas Delhomme & Thomas Riquelme & Iryna Shutava"
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
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#'#### Get the count data

#' # Raw data
#' ## Loading
#' Read the sample information
#' Create a table from the sequencing data (filename,sequencing date)
samples <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/rna_seq_spruce_needles_2011_2012_sample_info.csv")
w.samples <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/rna_seq_spruce_needles_winter_2012-2013_sample_info.csv")

#'List all the input kallisto tsv files and store them in orig variable
orig <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

winter <- list.files("../20170808/kallisto", 
             recursive = TRUE, 
             pattern = "abundance.tsv",
             full.names = TRUE)


#' name them
names(orig) <- sub("_sortmerna.*","",
                    sapply(strsplit(orig, "/"), .subset, 2))

names(winter) <- sub("_sortmerna.*","",
                     sapply(strsplit(winter, "/"), .subset, 4))

#' Combine the samples and the file list
samples <- rbind(samples,w.samples)
files <- c(orig,winter)

#' Match samples data to samples infos
files <- files[match(samples$Sequencing_ID,names(files))]
files <- files[!is.na(names(files))]

samples <- samples[samples$Sequencing_ID %in% names(files),]
stopifnot(all(samples$Sequencing_ID == names(files)))

#' The sample 1_120613_BC0TNNACXX_57_index18 has failed sequencing, so we
#' remove it
fail <- "1_120613_BC0TNNACXX_57_index18"
f.sel <- which(names(files)==fail)
samples <- samples[-f.sel,]
files <- files[-f.sel]

#' Extract the expression data from kallisto output files and store it into matrices
tx <- suppressMessages(tximport(files = files, 
                                type = "kallisto", 
                                txOut = TRUE))

# tximport(): looks inside the kallisto output files abundance.tsv, extract information and put it nicely into matrices for further use in downstream analysis
#store the matrices into tx
# str(tx)= a list of 4: abundance matrix, counts matrix, length matrix, countsFromAbundance (chr "no")

#' Store the counts in kt variable
kt <- round(tx$counts)
# round (keep only integer) and store the counts in kt 
# str(kt): a matrix
# rows: 66000 genes
# columns: 166 samples

#' Check for the genes that are never expressed
sel <- rowSums(kt) == 0 
# str(sel)= logical vector
# sums the counts for one gene in each samples
# show TRUE or FALSE for each gene if the gene is never expressed (has 0 counts in every time samples)

#' Check for the number of genes never expressed
sprintf("%s%% (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kt),digits=1),
        sum(sel),
        nrow(kt))
#in our experiment the number of genes and transcript are the same (one transcript for one gene)

#' Display the samples mean raw counts distribution
#' 
#' i.e. the mean raw count of every 
#' gene across samples is calculated
#' and displayed on a log10 scale.
#'
plot(density(log10(rowMeans(kt))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")
#density=likelihood, chance
# plot the likelihood that the mean of genes has this mean of raw counts 
# give an idea of the read depth, here the majority of gene maps ~30 reads (1.5 log10) 

#' Display all samples raw counts distribution
cols <- sample(pal,nlevels(samples$Sampling_Date),replace = TRUE)
plot.multidensity(mclapply(1:ncol(kt),function(k){log10(kt)[,k]},mc.cores=16L),
                  col=cols[as.integer(samples$Sampling_Date)],
                  legend.x="topright",
                  legend=levels(samples$Sampling_Date),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

# plot the likelihood that a gene in different time points has a certain number of counts 
# genes have a lot of chances (0.2 to 0.4 depending on the dates) to have a count of one 
# genes also have a lot of chances (0.2 to 0.3) to have between 2 and 100 counts 
# genes have fewer chances (0.2 to 0.0) to have between 100 and a 1000 counts per gene
# genes have even fewer chances (0.05 to 0.0) to have very high counts (1000 to 100000)= "powertail" 

#' The later samples (winter) have been sequence deeper than the older ones
dataset <- factor(ifelse(grepl("P7912",samples$Sequencing_ID),"latest","older"))
plot.multidensity(mclapply(1:ncol(kt),function(k){log10(kt)[,k]},mc.cores=16L),
                  col=pal[as.integer(dataset)],
                  legend.x="topright",
                  legend=levels(dataset),
                  legend.col=pal[1:2],
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE,recursive=TRUE)
write.csv(kt,file="analysis/kallisto/spruce-needles-unormalised-gene-expression_data.csv")
save(kt, samples, file = "analysis/kallisto/counts.rda")

#load("counts.rda")

#'#### Merge technical replicates

#' Merge technical replicates 

kt_tech_merged <- do.call(cbind,lapply(split.data.frame(t(kt),samples$Sampling_ID),colSums))
# transpose the count table kt to get the samples' names as rows in order to process them next (functions process rows and not columns)
# then, do the sums of counts for all genes at one date
# finally, bind again the splitted rows to reunified a count table with technical replicates merged for each dates 

#' Merge sampling dates for technical replicates
sampling_dates_merged <- factor(sapply(split(as.character(samples$Sampling_Date),
                                             samples$Sampling_ID),unique))

#' Save count table with technical replicates merged
save(sampling_dates_merged, kt_tech_merged, file="analysis/kallisto/kt_tech_merged.rda")

#' Plot
cols <- sample(pal,nlevels(sampling_dates_merged),replace = TRUE)
plot.multidensity(mclapply(1:ncol(kt_tech_merged),
                           function(k){log10(kt_tech_merged)[,k]},mc.cores=16L),
                  col=cols[as.integer(sampling_dates_merged)],
                  legend.x="topright",
                  legend=levels(sampling_dates_merged),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

#'#### Discard the NE samples which were resequenced with higher depth as TR00 
#'(keep TR00 with the same number)

#' find the bad samples
samples_names <- colnames(kt_tech_merged)
samples_names
new_names <- sub("NE-|TR00","",samples_names)
new_names
list <- split(samples_names, new_names)
list
duplicates <- list[sapply(list,length)==2]
bad_samples <- sapply(X = 1:length(duplicates),function(i){duplicates[[i]][1]})

#' remove from kt_tech_merged 
kt_tech_merged_cropped <- kt_tech_merged[,!colnames(kt_tech_merged) %in% bad_samples]

#' and remove the dates for this samples in the dates vector
sampling_dates_merged_cropped <- sampling_dates_merged[!names(sampling_dates_merged) %in% bad_samples]

# #' Save cropped dates and counts 
# Save(kt_tech_merged_cropped, sampling_dates_merged_cropped, file="kt_tech_merged_cropped.rda")

#' order dates chronologically
dates <- sampling_dates_merged_cropped
ordered_dates <- dates[order(as.Date(dates, format="%Y-%m-%d"))]

#' order count table chronologically
kt_tech_merged_cropped_ordered <- kt_tech_merged_cropped[,names(ordered_dates)]

#' Save ordered dates and counts
save(ordered_dates, kt_tech_merged_cropped_ordered, file = "kt_tech_merged_cropped_ordered.rda")

plot.multidensity(mclapply(1:ncol(kt_tech_merged_cropped_ordered),
                           function(k){log10(kt_tech_merged_cropped_ordered)[,k]},mc.cores=16L),
                  col=cols[as.integer(ordered_dates)],
                  legend.x="topright",
                  legend=levels(ordered_dates),
                  legend.col=cols,
                  legend.lwd=2,
                  main="samples raw counts distribution",
                  xlab="per gene raw counts (log10)")

#'##### Data normalisation 

#' ### vst blind

#load technical replicates merged data cropped for NE bad samples
load("kt_tech_merged_cropped_ordered.rda")

#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate
#'
#' Create the dds object, without giving any prior on the design

dds.kt <- DESeqDataSetFromMatrix(
  countData = kt_tech_merged_cropped_ordered,
  colData = data.frame(date=ordered_dates),
  design = ~ date)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kt <- estimateSizeFactors(dds.kt)
sizes.kt <- sizeFactors(dds.kt)
# names(sizes.kt) <- colnames(kt)
# not useful sizes.kt already have the right names
pander(sizes.kt)
boxplot(sizes.kt, main="Sequencing libraries size factor")
abline(h=1,lty=2,col="gray")

#' ## Variance Stabilising Transformation
#' ### Blind
system.time(vsd.kt <- varianceStabilizingTransformation(dds.kt, blind=TRUE))
vst_blind <- assay(vsd.kt)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_blind <- vst_blind - min(vst_blind) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_blind[rowSums(vst_blind)>0,]) #mean variance stabilized between 0.5 and 1
meanSdPlot(log2(kt+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kt,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")
#VST validated: mean variance stabilized around 0.5

write.csv(vst_blind,"analysis/kallisto/spruce-needles-normalized-data-vst-blind.csv")

save(ordered_dates, vst_blind, file="vst_blind.rda")

#' ### vst aware (calculated accordingly to samplingDates)

#' Calculate vst aware with DESeq2 taking the date as a parameter
#' ### Model-aware
system.time(vsd2 <- varianceStabilizingTransformation(dds.kt, blind=FALSE))
vst_aware <- assay(vsd2)

#colnames(vst.kt) <- colnames(kt)
#already have the right names
vst_aware <- vst_aware - min(vst_aware) 
# substracts from all values the smaller value, so the smaller value is now 0

#' Validate the VST 
meanSdPlot(vst_aware[rowSums(vst_aware)>0,]) #mean variance stabilized between 0 and 0.5
meanSdPlot(log2(kt+1)) #without Variance Stabilizing Transformation: variation around 1
meanSdPlot(log2(counts(dds.kt,normalized=TRUE)+1)) #normalized read depth: variance is not stabilized ("curvy")

# vst_aware validated: mean variance is way lower when vst is calculated "aware of the date factor" compare to before when vst was calculated "blind"
# stabilized around 0.2 instead of 0.5 when vst was blind

write.csv(vst_aware,"analysis/kallisto/spruce-needles-normalised-data-vst-aware.csv")
save(ordered_dates, vst_aware, file="vst_aware.rda")
load("vst_aware.rda")

#' ### Filter data for expressed enough genes 

## REPLACE the function by a call to the featureSelection.R utility - check in arabidopsis pib project + KallistoBiologicalQA.R (from line 197) 

## IS: ----------------------

#' ## Select the genes that are expressed 
# vst_aware

sels <- sapply(1:10,function(i){
  featureSelect(vst_aware,conditions = factor(paste0(samples$Sampling_ID,samples$Sampling_Date)),
                exp=i)})


plot(colSums(sels),type="l",xlab="vst cutoff",
     main="number of genes selected at cutoff",ylab="number of genes")

sel <- sels[,2]



#' Hierarchical clustering of the data
#plot(hclust(dist(t(vst_aware[sel,]))),labels=samples$Sampling_ID)


#' Create a heatmap
#hpal <- colorRampPalette(c("blue","white","red"))(100)

#heatmap.2(vst_aware[sel,],trace="none",col=hpal)

#s.vst <- t(scale(t(vst_aware)))

#library(hyperSpec)

#heatmap.2(s.vst[sel,],distfun = pearson.dist,
#          hclustfun = function(X){hclust(X,method="ward.D")},
#          trace="none",col=hpal,labRow = FALSE,
#          labCol = paste(samples$Genotype,samples$Time,sep="-"))


#hc <- hclust(pearson.dist(s.vst[sel,]),method = "ward.D")

#tc <- cutree(hc,k=4)

#tc

#nams <- split(names(tc),tc) 
#library(IRanges)
#barplot(elementNROWS(nams))

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
