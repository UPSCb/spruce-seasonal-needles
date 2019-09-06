#' ---
#' title: "Spruce Needles Kallisto VST model-aware"
#' author: "Nicolas Delhomme & Thomas Riquelme"
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

#' # Data normalisation 

#' load technical replicates merged data cropped for NE bad samples
load("analysis/kallisto/kt_tech_merged_cropped_ordered.rda")

#' Create the dds object, without giving any prior on the design
dds.kt <- DESeqDataSetFromMatrix(
  countData = kt_tech_merged_cropped_ordered,
  colData = data.frame(date=ordered_dates),
  design = ~ date)

#' ## VST aware (calculated accordingly to samplingDates)
vsd2 <- system.time(varianceStabilizingTransformation(dds.kt, blind=FALSE))

#' # Export
save(ordered_dates, vsd2, file="analysis/kallisto/vst_aware.rda")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
