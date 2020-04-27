#' ---
#' title: "Plot Expression Profile"
#' author: "Thomas Riquelme, Iryna Shutava & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Libraries
suppressPackageStartupMessages({
  library(amap)
  library(here)
  library(corrplot)
  library(ggplot2)
  library(matrixStats)
})

#' * Data
#' 
#' table of counts vst transformed
load(here("data/vst_aware.rda"))

#' Focus on one year of sampling - samples past 2012-04-24 are the
#' next generation of needles
sel <- ordered_dates %in% levels(ordered_dates)[1:(grep("2012-05",levels(ordered_dates))[1]-1)]
ordered_dates <- factor(as.character(ordered_dates)[sel])
vst_aware <- vst_aware[,sel]

#' Change the dates nomenclature
ddf <- do.call(rbind,strsplit(levels(ordered_dates),"-"))
labels <- ordered_dates
levels(labels) <- ifelse(as.integer(ddf[,3]) >= 25,
                         month.name[as.integer(ddf[,2])+1],
                         ifelse(as.integer(ddf[,3]) <= 10,
                                month.name[as.integer(ddf[,2])],
                                paste0("mid-",month.name[as.integer(ddf[,2])])
                         ))

plotLabels <- as.character(labels[match(levels(ordered_dates),ordered_dates)])
plotLabels[duplicated(plotLabels)] <- ""

#' * Combine the replicates
vst_bio_rep_mean <- do.call(cbind,lapply(split.data.frame(t(vst_aware),ordered_dates),colMeans))

vst_bio_rep_median <- do.call(cbind,lapply(split.data.frame(t(vst_aware),ordered_dates),colMedians))
rownames(vst_bio_rep_median) <- rownames(vst_bio_rep_mean)

vst_bio_rep_sd <- do.call(cbind,lapply(split.data.frame(t(vst_aware),ordered_dates),colSds))

vst_bio_rep_lwr <- vst_bio_rep_mean-qt(0.975,3)*vst_bio_rep_sd/sqrt(3)

#' "2011-06-30" has no replicates, thus it is impossible to calculate sd, lwr and upr
#' replacement of NA values  in upr and lwr by the mean to be able to plot the confidence interval of the mean later
vst_bio_rep_lwr[,"2011-06-30"] <- vst_bio_rep_mean[,"2011-06-30"]
vst_bio_rep_lwr[vst_bio_rep_lwr < 0] <- 0

vst_bio_rep_upr <- vst_bio_rep_mean+qt(0.975,3)*vst_bio_rep_sd/sqrt(3)

vst_bio_rep_upr[,"2011-06-30"] <- vst_bio_rep_mean[,"2011-06-30"]

#' Scale the data
s.vst <- t(scale(t(vst_aware)))

s.vst.mean <- do.call(cbind,lapply(split.data.frame(t(s.vst),ordered_dates),colMeans))

s.vst.median <- do.call(cbind,lapply(split.data.frame(t(s.vst),ordered_dates),colMedians))
rownames(s.vst.median) <- rownames(s.vst.mean)

s.vst.sd <- do.call(cbind,lapply(split.data.frame(t(s.vst),ordered_dates),colSds))

s.vst.lwr <- s.vst.mean-qt(0.975,3)*s.vst.sd/sqrt(3)

#' "2011-06-30" has no replicates, thus it is impossible to calculate sd, lwr and upr
#' replacement of NA values  in upr and lwr by the mean to be able to plot the confidence interval of the mean later
s.vst.lwr[,"2011-06-30"] <- s.vst.mean[,"2011-06-30"]

s.vst.upr <- s.vst.mean+qt(0.975,3)*s.vst.sd/sqrt(3)

s.vst.upr[,"2011-06-30"] <- s.vst.upr[,"2011-06-30"]

#' * Gene Of Interest
goi <- read.csv(here("doc/GOI_list.csv"), header = FALSE,as.is=TRUE)

#' read GOI names file
goi_names <- read.csv(here("doc/GOI_names.txt"), header=FALSE)

#' 2020-03-25: Add new extra gois
extra.goi <- data.frame(c("MA_734271g0010.1", "MA_210976g0010.1", 
               "MA_60075g0010.1", "MA_90810g0010.1"))
goi <- rbind(goi,extra.goi)
goi_names <- c(goi_names,extra.goi)


#' ```{r test plot, eval=FALSE, echo=FALSE}
#' # Test to plot one gene of the list
#' psbs_mean <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) == goi[52,1]]
#' psbs_median <- vst_bio_rep_median[rownames(vst_bio_rep_median) == goi[1,1]]
#' psbs_lwr <- vst_bio_rep_lwr[rownames(vst_bio_rep_lwr) == goi[52,1]]
#' psbs_upr <- vst_bio_rep_upr[rownames(vst_bio_rep_upr) == goi[52,1]]
#' 
#' DF <- data.frame(time=as.Date(levels(ordered_dates)),
#'                  mean=psbs_mean,
#'                  lwr=psbs_lwr,
#'                  upr=psbs_upr)
#' 
#' p <- ggplot(DF, aes(time, group = 1)) +
#'   geom_line(aes(y=mean), color="blue") +
#'   geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) +
#'   scale_x_date(date_breaks = "month") +
#'   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#'   labs(y = "vst counts", title = "New plot title") +
#'   ggtitle("Expression of gene test", subtitle = NULL)
#' plot(p)
#' 
#' 
#' psbs_mean <- s.vst.mean[rownames(s.vst.mean) == goi[52,1]]
#' psbs_median <- s.vst.median[rownames(s.vst.median) == goi[1,1]]
#' psbs_lwr <- s.vst.lwr[rownames(s.vst.lwr) == goi[52,1]]
#' psbs_upr <- s.vst.upr[rownames(s.vst.upr) == goi[52,1]]
#' 
#' DF <- data.frame(time=as.Date(levels(ordered_dates)),
#'                  mean=psbs_mean,
#'                  lwr=psbs_lwr,
#'                  upr=psbs_upr)
#' 
#' p <- ggplot(DF, aes(time, group = 1)) +
#'   geom_line(aes(y=mean), color="blue") +
#'   geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) +
#'   scale_x_date(date_breaks = "month") +
#'   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#'   labs(y = "standard score", title = "New plot title") +
#'   ggtitle("Expression of gene test", 
#'           subtitle = paste("Average vst expression:",
#'                            round(mean(vst_bio_rep_mean[rownames(vst_bio_rep_mean) == goi[52,1]]),digits=2)))
#' plot(p)
#' 
#' 
#' 
#' ```

#' Do a loop to plot the expression profile for each gene of the list
for (i in 1:nrow(goi)) {
  message(i)
  #check if goi are in vst with %in%
  if (goi[i,1] %in% rownames(vst_aware)) {
    g_median <- vst_bio_rep_median[rownames(vst_bio_rep_median) == goi[i,1]]
    g_mean <- vst_bio_rep_mean[rownames(vst_bio_rep_mean) == goi[i,1]]
    g_lwr <- vst_bio_rep_lwr[rownames(vst_bio_rep_lwr) == goi[i,1]]
    g_upr <- vst_bio_rep_upr[rownames(vst_bio_rep_upr) == goi[i,1]]
    DF <- data.frame(time=as.Date(levels(ordered_dates)),
                     mean=g_median,
                     lwr=g_lwr,
                     upr=g_upr)
    
    p <- ggplot(DF, aes(time, group = 1)) +
      annotate("rect", xmin=as.Date("2011-10-01"),
                    xmax=as.Date("2012-03-01"),
                    ymin=-Inf,ymax=+Inf,
                    fill="#ADD8E6",alpha=0.3) +
      geom_line(aes(y=mean), color="blue") +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) + 
      scale_x_date(date_breaks = "month") +
      theme_bw() + coord_cartesian(ylim=c(0,max(vst_bio_rep_upr))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      labs(y = "vst counts", title = paste("Expression profile of", goi_names[i,1]))
    plot(p)
    suppressMessages(ggsave(filename = here("data/expression_profiles",
                                            paste0(goi_names[i,1],"_expression_profile.jpeg"))))
    
    DF <- data.frame(time=as.Date(levels(ordered_dates)),
                     mean=s.vst.mean[rownames(s.vst.mean) == goi[i,1]],
                     lwr=s.vst.lwr[rownames(s.vst.lwr) == goi[i,1]],
                     upr=s.vst.upr[rownames(s.vst.upr) == goi[i,1]])
    
    p <- ggplot(DF, aes(time, group = 1)) +
      annotate("rect", xmin=as.Date("2011-10-01"),
               xmax=as.Date("2012-03-01"),
               ymin=-Inf,ymax=+Inf,
               fill="#ADD8E6",alpha=0.3) +
      geom_line(aes(y=mean), color="blue") +
      geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.3) +
      scale_x_date(date_breaks = "month") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(y = "standard score") +
      ggtitle(paste("Scaled expression profile of", goi_names[i,1]), 
              subtitle = paste("Average vst expression:",
                               round(mean(g_mean),digits=2)))
    plot(p)
    suppressMessages(ggsave(filename = here("data/expression_profiles",
                                            paste0(goi_names[i,1],"_scaled_expression_profile.jpeg"))))
    
  }
}

#' Remark ELIP_C does not match our data genes id

#'### Hierarchical clustering to see which genes have similar patterns

#' Improve gene of interest data (add ELIP_C)
goi <- cbind(goi,goi_names)
colnames(goi) <- c("id","name")

#' select mean vst counts (of biological replicates for one time point) for gene of interest (goi)
sel <- match(goi$id, rownames(vst_bio_rep_mean))
vst_bio_rep_mean_goi <- vst_bio_rep_mean[sel,]

#' replace gene ids by gene names as rownames
rownames(vst_bio_rep_mean_goi) <- goi$name

#' remove ELIP_C
goi <- goi[! is.na(sel),]
vst_bio_rep_mean_goi <- vst_bio_rep_mean_goi[! is.na(sel),]

#' Normalize data to obtain z-score to quantify only the variation around the mean (=pattern) and not the amplitude anymore for each each gene
#' by doing so we can compare the genes relatively to their expression pattern and not relatively of their amplitude like before
vst_bio_rep_mean_goi_scaled <- t(scale(t(vst_bio_rep_mean_goi)))

#' Remove RabA1, which is not expressed
sel <- which(rowSums(is.na((vst_bio_rep_mean_goi_scaled))) > 0)
vst_bio_rep_mean_goi_scaled <- vst_bio_rep_mean_goi_scaled[-sel,]
vst_bio_rep_mean_goi <- vst_bio_rep_mean_goi[-sel,]

#' qqplot to check if our genes expression follows a Normal distribution or not
qqnorm(vst_bio_rep_mean_goi)
qqline(vst_bio_rep_mean_goi, col=3)
qqnorm(vst_bio_rep_mean_goi_scaled)
qqline(vst_bio_rep_mean_goi_scaled, col = 2)
#' ==> does not follow a line, the distribution is not Normal
#' ==> then Spearman correlation should be use because is not parametric

#' Perform hierarchical clustering
#' Firstly compute distance according to "spearman correlation" method
#' Secondly compute clustering according to "complete" linkage
vst_bio_rep_mean_goi_cluster <- hcluster(vst_bio_rep_mean_goi, method = "spearman", link = "complete")
vst_bio_rep_mean_goi_scaled_cluster <- hcluster(vst_bio_rep_mean_goi_scaled, method = "spearman", link = "complete")

#' Plot the dendrogram
plot(vst_bio_rep_mean_goi_cluster, 
     main="Cluster Dendrogram of photosynthetic genes",
     cex=0.8)
plot(vst_bio_rep_mean_goi_scaled_cluster, 
     main="Cluster Dendrogram of photosynthetic genes",
     cex=0.8)

pdf(file=here("data/photosyntetic-genes-dendrogram.pdf"),width=12,height=8)
plot(vst_bio_rep_mean_goi_scaled_cluster, 
     main="Cluster Dendrogram of photosynthetic genes",
     cex=0.8,xlab="photosyntetic genes",sub="spearman distance, complete linkage")
dev.off()

#' Plot the correlation
pdf(file = here("data/expression_profiles/photosynthetic_genes_correlation_matrix.pdf"))
corrplot(cor(t(vst_bio_rep_mean_goi)),method="square",order="hclust", 
         title = "Correlation matrix of photosynthetic genes",
         tl.cex=0.4, cl.cex=0.5, tl.col="black", addrect=2, is.corr = FALSE, mar = c(0,0,2,0))
dev.off()
corrplot(cor(t(vst_bio_rep_mean_goi_scaled)),method="square",order="hclust",
        tl.cex=0.4, cl.cex=0.5, tl.col="black", addrect=2, is.corr = FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

