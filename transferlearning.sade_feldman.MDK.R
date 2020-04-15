#############################################################################
# Transfer Learning Analysis: Use projectR package from Fertig Lab to project
# CoGAPs laten factors (i.e. patterns) from a single-cell mouse
# immune-treatment dataset to a single-cell human immune treatment dataset
# In other words, learn the representation of each pattern learned by CoGAPs
# in the mouse data set within the human dataset by using projectR
# Author: Emily Davis Marcisak
# Updated by Michael D. Kessler
#############################################################################

#############################################################################
# Set up environment
#############################################################################

# load packages
library(projectR)
library(CePa)
library(org.Hs.eg.db)
library(biomaRt)
library(gplots)
library(reshape2)
library(ggplot2)
library(CoGAPS)
library(data.table)
library(ComplexHeatmap)
library(viridis)
library(GEOquery)
library(RColorBrewer)
library(ROCR)
library(dplyr)
library(ggalluvial)
library(plyr)

# set working directory
setwd('/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI')

#############################################################################
# Download data - you only need to download data once
#############################################################################

# download data from GEO using GEOquery function
#filePaths = getGEOSuppFiles("GSE120575")

#############################################################################
# Loading and preprocessing datasets
#############################################################################

## load human single cell dataset (sade-feldman) - TPM data
# Note: these files are quite big (~4GB). Even though they can be kept to
# ~200MB using gzip compression (probably due to sparsity), they are expanded
# very considerably and may cause memory trouble in R/RStudio. If using
# RStudio in particular, make sure to update your RAM limits according to the
# following suggestion https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
# This should prevent memory errors and faciliate a smooth and reproducible
# workflow

# if you downloaded data with command above with this session,
# you can use the following (currently commented out) command to load data.
# Otherwise, load using an explicit path, as seen below.
#sfcounts <- data.table::fread(file = rownames(filePaths)[1], fill = TRUE)
sfcounts <- data.table::fread(file = "GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz", fill = TRUE)
# dim(sfcounts) # check dimensions, which should be 55739 16293

## save barcodes, which will be necessary to link to metadata
barcodes <- sfcounts[c(1),]
## remove extra metadata rows
sfcounts <- sfcounts[-c(1:2),]
## check out format again, looking at end
#View(sfcounts[1:5,16289:16293])
## remove extra column at end, which seems to be filled with NAs 
sfcounts <- sfcounts[,-c(16293)]
sfgenes <- sfcounts[, 1]$V1 # save gene symbols for reassignment below
# set as matrix, as current data type is strange and naming rows fails
sfcounts <- as.matrix(sfcounts)
## now remove extra name column
sfcounts <- sfcounts[,-c(1)]
# convert to numeric character type
sfcounts <- as.matrix(apply(sfcounts, 2, as.numeric))
# set rownames
rownames(sfcounts) <- sfgenes

## load sade-feldman meta data
sfmeta <- read.table(file = "GSE120575/GSE120575_patient_ID_single_cells.txt.gz", sep = '\t', header = FALSE)
# proess 
## check end
#View(sfmeta[16306:16343,])
## cut out extra rows at end
sfmeta <- sfmeta[c(1:16306),]
## cut out extra cols at end
sfmeta <- sfmeta[,-c(12:35)]
## cut out extra beginning rows
sfmeta <- sfmeta[-c(1:15),]
## check
#View(sfmeta[1:5,1:5])
# manually set colnames
colnames(sfmeta) <- c("Sample name", "title", "source name", "organism", "patient ID (Pre=baseline; Post= on treatment)", 
                     "response", "therapy")
## remove remaining empty cols
sfmeta <- sfmeta[,-c(8:14)]
## set sfcounts colnames to barcode names from meta data
colnames(sfcounts) <- sfmeta$title
## split and combine meta information for grouping
split <- strsplit(as.character(sfmeta[,5]),split='_', fixed=TRUE)
p1 <- sapply(split, "[", 1)
p2 <- sapply(split, "[", 2)
sfmeta$patient <- p2
sfmeta$status <- p1
sfmeta$combined <- paste(sfmeta$status, sfmeta$therapy, sfmeta$response, sep = " ")
## convert gene names to the same format as the CoGAPS object
symbols <- toupper(rownames(sfcounts))

## checking summary table of treatments
#table(sfmeta$combined)

#############################################################################
# Run transfer learning
###################################################
## load CoGAPS result object
load("myeloid_run7_result3.RData")
gaps2mel <- projectR(data = sfcounts, 
                      loadings = gapsResult,
                      full = TRUE, 
                      dataNames = symbols)

## get projected matrix out of new TL object
TL.proj <- gaps2mel[["projection"]]

## add back in sample or treatment names
colnames(TL.proj) <- sfmeta$combined
sfmeta <- droplevels(sfmeta) # drop unused factor levels
## set any NAs to zero
TL.proj[which(is.na(TL.proj))] <- 0
## confirm NAs were removed and no other atypical value types (Inf, Nan, etc)
# which(is.na(TL.proj))
# which(is.infinite(TL.proj))
# which(is.nan(TL.proj))

## save TL results
#saveRDS(TL.proj, file = "human_ici_scrna_projection.rds")
# read in if starting analysis from here and don't want to rerun steps above
TL.proj <- readRDS("human_ici_scrna_projection.rds")

############################################################################
# Graph results
#############################################################################

# Plot Complex Heatmaps
# simplify response coding for use as row annotation - use colnames for this
simple.row_anno <- gsub("Pre ", "", colnames(TL.proj))
simple.row_anno <- gsub("Post ", "", simple.row_anno)
simple.row_anno <- as.factor(simple.row_anno)
# set cols - use viridis palette
#mycols <- viridis(length(levels(simple.row_anno)), option = "D")
mycols <- rainbow(6)
# open connection to plotting device
tiff(file = "TL.all_CoGAPs_patterns.ComplexHeatmap.pre_post.tiff",
     width = 6000, height = 5000, units = "px", res = 800)
row_ha = rowAnnotation(Response = simple.row_anno,
            col = list(Response = c(
            "anti-CTLA4 Non-responder" = mycols[1],
            "anti-CTLA4 Responder" = mycols[2],
            "anti-CTLA4+PD1 Non-responder" = mycols[3],
            "anti-CTLA4+PD1 Responder" = mycols[4],
            "anti-PD1 Non-responder" = mycols[5],
            "anti-PD1 Responder" = mycols[6]

)))
# scale by sample
TL.proj.rownames <- rownames(TL.proj)
TL.proj.scl <- apply(TL.proj, 2, scale)
rownames(TL.proj.scl) <- TL.proj.rownames # reassign col as row names
# call heatmap from complexheatmap package
ComplexHeatmap::Heatmap(t(TL.proj.scl), col=inferno(100), name = "mat",
        clustering_distance_rows = "pearson",
        column_names_gp = gpar(fontsize = 5),
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 2,
        show_row_names = FALSE
)
dev.off()

# now split up the pre and post treatment samples and see if that improves
# the visualization of clustering
TL.proj.pre <- TL.proj[, which(sapply(colnames(TL.proj), function(x) grepl("Pre", x)))]
TL.proj.post <- TL.proj[, which(sapply(colnames(TL.proj), function(x) grepl("Post", x)))]

# plot pre samples
# open connection to plotting device
tiff(file = "TL.all_CoGAPs_patterns.ComplexHeatmap.pre.tiff",
     width = 6000, height = 5000, units = "px", res = 800)
row_ha = rowAnnotation(Response = colnames(TL.proj.pre),
                       col = list(Response = c(
                         "Pre anti-CTLA4 Non-responder" = mycols[1],
                         "Pre anti-CTLA4 Responder" = mycols[2],
                         "Pre anti-CTLA4+PD1 Non-responder" = mycols[3],
                         "Pre anti-CTLA4+PD1 Responder" = mycols[4],
                         "Pre anti-PD1 Non-responder" = mycols[5],
                         "Pre anti-PD1 Responder" = mycols[6]
                       )))
# scale by sample
TL.proj.pre.rownames <- rownames(TL.proj.pre)
TL.proj.pre.scl <- apply(TL.proj.pre, 2, scale)
rownames(TL.proj.pre.scl) <- TL.proj.pre.rownames # reassign col as row names
# call heatmap from complex heatmap package
ComplexHeatmap::Heatmap(t(TL.proj.pre.scl), col=inferno(100), name = "mat",
        clustering_distance_rows = "pearson",
        column_names_gp = gpar(fontsize = 5),
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 2,
        show_row_names = FALSE
)
dev.off()

# plot post samples
# open connection to plotting device
tiff(file = "TL.all_CoGAPs_patterns.ComplexHeatmap.post.tiff",
     width = 6000, height = 5000, units = "px", res = 800)
row_ha = rowAnnotation(Response = colnames(TL.proj.post),
                       col = list(Response = c(
                         "Post anti-CTLA4+PD1 Non-responder" = mycols[6],
                         "Post anti-CTLA4+PD1 Responder" = mycols[1],
                         "Post anti-PD1 Non-responder" = mycols[4],
                         "Post anti-PD1 Responder" = mycols[5]
                       )))
# scale by sample
TL.proj.post.rownames <- rownames(TL.proj.post)
TL.proj.post.scl <- apply(TL.proj.post, 2, scale)
rownames(TL.proj.post.scl) <- TL.proj.post.rownames # reassign col as row names
# call heatmap from complex heatmap package
ComplexHeatmap::Heatmap(t(TL.proj.post.scl), col=inferno(100), name = "mat",
        clustering_distance_rows = "pearson",
        column_names_gp = gpar(fontsize = 5),
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 2,
        show_row_names = FALSE
)
dev.off()

# Now plot pattern distributions by treatment_response groups
# using violin plots and boxplots

# First do some preprocessing
# read in single cell annotations
sfanno <- read.csv("sfanno.csv", header = T)
# transpose TL projection data and convert to df
TL.proj.df = as.data.frame(t(TL.proj))
# add single cell annotatinos to TL.proj.df
mtch <- match(sfmeta$title, sfanno$Cell.Name)
TL.proj.df$cluster <- as.character(sfanno[mtch,]$Cluster.number)
# Also set cluster factors to cell types
TL.proj.df$cluster <- plyr::revalue(TL.proj.df$cluster, c("1" = "B-cells",
                                    "2" = "Plasma cells",
                                    "3" = "Monocytes/Macrophages",
                                    "4" = "Dendritic cells",
                                    "5" = "Lymphocytes",
                                    "6" = "Exhausted CD8+ T-cells",
                                    "7" = "Regulatory T-cells",
                                    "8" = "Cytotoxicity",
                                    "9" = "Exhausted/HS CD8+ T-cells",
                                    "10" = "Memory T-cells",
                                    "11" = "Lymphocytes exhausted/cell-cycle"))
TL.proj.df$group <- rownames(TL.proj.df)
TL.proj.df$status <- sfmeta$status
TL.proj.df.m <- reshape2::melt(TL.proj.df, varnames = c("cluster", "group",
      "status"), na.rm = FALSE, as.is = FALSE, value.name = "value")
# reformat group labels
colnames(TL.proj.df.m) <- c("cluster", "group", "status", "pattern", "weight")
TL.proj.df.m$group <- gsub(".", " ", TL.proj.df.m$group, fixed=T)
TL.proj.df.m$group <- gsub("[[:digit:]]", "", TL.proj.df.m$group)
TL.proj.df.m$group <- gsub("CTLA", "CTLA4", TL.proj.df.m$group)
TL.proj.df.m$group <- gsub("ponder ", "ponder", TL.proj.df.m$group)
TL.proj.df.m$group <- as.factor(TL.proj.df.m$group)

## violin plot of projected pattern weights
num_pats <- rownames(TL.proj)
plot_list = list()
for(i in 1:length(num_pats)){
  pltdat <- TL.proj.df.m[TL.proj.df.m$pattern == num_pats[i],]
  plotname <- paste("Pattern", i, "of", num_pats[i], sep = " ")
  p <- ggplot(pltdat, aes(x=pattern, y=weight, fill = group)) +
    geom_violin(alpha=0.8) + theme_bw() +
    theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
          axis.ticks.x=element_blank()) +
    scale_fill_viridis_d() + ggtitle(plotname) +
    geom_point(alpha = 0.8, aes(fill = group), size = 1.5, shape = 21, position = position_jitterdodge()) +
    ylab("pattern weight")
  plot_list[[i]] = p
}
# save plots together, across PDF pages
pdf("21pats_TL_sc_og_violinplots_trim.pdf", width = 7, height = 6)
for (i in 1:length(num_pats)) {
  print(plot_list[[i]])
}
dev.off()

# make faceted boxplots
pdf("TL.all_CoGAPs_patterns.faceted_boxplots.pdf", width = 18, height = 14)
ggplot(TL.proj.df.m, aes(x=group, y=weight, fill = group)) +
  geom_boxplot(alpha=0.8, notch = T) + theme_bw() +
  theme(axis.title=element_text(size=30),
        axis.text.x=element_blank(), #element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 20)) +
  scale_fill_viridis_d() + ggtitle(NULL) +
  facet_wrap(~pattern, scales = "free") +
  theme(legend.position="bottom",
        legend.text=element_text(size=15))
  #geom_point(alpha = 0.8, aes(fill = group), size = 1.5, shape = 21, position = position_jitterdodge()) +
  #ylab("pattern weight")
dev.off()

# remake notched boxplots after splitting cells by single cell annotations
# per pattern, make faceted bocplots per single cell annotation, with each
# panel having all treatment response groups
for (i in 1:21){
  # pp stands for per pattern
  TL.proj.df.m.pp <- subset(TL.proj.df.m, pattern == paste0("Pattern_", i))
  pdf(paste0("dists.per_pattern.per_cell/TL.pattern", i, ".faceted_violinplots.pdf"), width = 18, height = 14)
  plt <- ggplot(TL.proj.df.m.pp, aes(x=group, y=weight, fill = group)) +
    geom_violin(alpha=0.8) + theme_bw() +
    theme(axis.title=element_text(size=30),
          axis.text.x=element_blank(), #element_text(size=rel(0.75), angle=45, hjust = 1),
          axis.ticks.x=element_blank(),
          strip.text.x = element_text(size = 16)) +
    scale_fill_viridis_d() + ggtitle(NULL) +
    facet_wrap(~cluster, scales = "free") +
    theme(legend.position="bottom",
          legend.text=element_text(size=15))
  #geom_point(alpha = 0.8, aes(fill = group), size = 1.5, shape = 21, position = position_jitterdodge()) +
  #ylab("pattern weight")
  print(plt)
  dev.off()
}
