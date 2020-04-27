#' ---
#' title: "Plot Expression Profile"
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
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(plotly)
})
  
#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#' Source multidensity function
#source("~/Git/UPSCb/src/R/plot.multidensity.R")

load(file = "vst_aware.rda")

#' remove the duplicated dates to use as plot legend 
dates <- ordered_dates
dates[duplicated(dates)] <- ""
factor_dates <- factor(dates)
factor_dates

#' Store data about gene of interest psbs relative to photosynthesis
psbs1 <- vst_aware[rownames(vst_aware) == "MA_575444g0010.1"]
psbs2 <- vst_aware[rownames(vst_aware) == "MA_571554g0010.1"]
psbs3 <- vst_aware[rownames(vst_aware) == "MA_4107548g0010.1"]

#' Plot the expression of all 3 psbs genes

plotmatrix(rbind(psbs1,psbs2,psbs3), 
           ylab = "counts (vst)", 
           las=2,
           xlabels = factor_dates, 
           cols = c("blue", "green", "red"),
           main = "expression profiles of psbs genes")

#read Gene Of Interest csv file
goi <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_list.csv", header = FALSE)

#write.csv(x = goodid, file = "~/Git/UPSCb/projects/spruce-needles/doc/goodid.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
#read GOI names file
goi_names <- read.csv("~/Git/UPSCb/projects/spruce-needles/doc/GOI_names.txt", header=FALSE)

# Test to plot one gene of the list
j <- vst_aware[rownames(vst_aware) == goi[1,1]]

# jpeg(filename = "Rplot%03d.jpeg",
# width = 480, height = 480, units = "px", pointsize = 12,
# quality = 75,
# bg = "white", res = NA, ...,

jpeg(filename = "test.jpg",
     width = 1200, height = 736)
par(mar=c(6.1,4.1,4.1,2.1))
plot(j, 
     type = "l",
     las = 2,
     ylab = "counts (vst)",
     xaxt = "n",
     xlab = NA,
     col ="blue",
     main = "expression profile of psbs"
)
axis(1, at=1:90, labels=factor_dates, las = 2)

#dev.print(jpeg,'psbs3.jpeg')
dev.off()

par(mar=mar)

#plot(psbs2)
#plot(psbs3,type="l")

# do a loop to plot an expression profile for each gene of the list
for (i in 1:nrow(goi)) {
  message(i)
  #check if goi are in vst with %in%
  if (goi[i,1] %in% rownames(vst_aware)) {
    jpeg(filename = paste0(goi_names[i,1],"_expression_profile.jpeg"), width = 1200, height = 736)
    plot(vst_aware[rownames(vst_aware) == goi[i,1]], 
         type = "l",
         las = 2,
         ylab = "counts (vst)",
         xaxt = "n",
         xlab = NA,
         col ="blue",
         main = paste("Expression profile of", goi_names[i,1])
    )
    axis(1, at=1:90, labels=factor_dates, las = 2)
    dev.off()
  }
}

# genes CytC6 and FC1 did not match with vst genes because they were not in Congenie database
# "only found with BLAST, but not through Congenie search function"
# ELIP_C does not match either our data genes id
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
