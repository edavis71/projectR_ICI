#############################################################################
# Transfer Learning Analysis: Use projectR package from Fertig Lab to project
# CoGAPs laten factors (i.e. patterns) from a single-cell mouse
# immune-treatment dataset to TCGA data for each cancer type
# In other words, learn the representation of each pattern learned by CoGAPs
# in the mouse data set within TCGA cancer data and look for correlations
# between patterns and survival
# Authors: Michael D. Kessler and Emily Davis Marcisak
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
library(TCGAbiolinks)
library(SummarizedExperiment)
# set working directory
setwd('/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI')

###################################################
### download data per cancer type
###################################################

# vector of cancer type names
cancer_types <- c("LAML", "ACC","BLCA", "LGG", "BRCA", "CESC", "CHOL",
                  "LCML", "COAD", "CNTL", "ESCA", "FPPP", "GBM", "HNSC",
                  "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "DLBC",
                  "MESO", "MISC", "OV", "PAAD", "PCPG", "PRAD", "READ",
                  "SARC", "SKCM", "STAD", "TGCT","THYM", "THCA", "UCS",
                  "UCEC")
# format query labels
TCGA_cancer_labels <- paste0("TCGA-", cancer_types)
# download data per cancer type using TCGABioLinks
TCGA_data <- list() # list of tissue specific TCGA data frames
for (ct in TCGA_cancer_labels){
  # use TCGABioLinks to download data and parse into dataframes
  # define data you want to download
  query <- GDCquery(project = ct,
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  # download data
  GDCdownload(query, method = "api", files.per.chunk = 10)
  # parse/prepare the data 
  TCGAdata <- GDCprepare(query)
  # store cancer type specific meta data
  metaTCGA <- as.data.frame(SummarizedExperiment::colData(TCGAdata))
  # store cancer type specific expression data
  expTCGA <- as.matrix(SummarizedExperiment::assay(TCGAdata))
  # store this data in the TCGA_data list
  TCGA_data[[ct]] <- list()
  TCGA_data[[ct]][["meta"]] <- metaTCGA
  TCGA_data[[ct]][["exp"]] <- expTCGA
}

#############################################################################
# combine cancer type specific meta and expression data frames
#############################################################################

metaTCGA <- NULL # reinit
expTCGA <- NULL # reinit
for (ct in TCGA_cancer_labels){
  # extract cancer speicfic metadata
  tempmeta <- TCGA_data[[ct]][["meta"]]
  tempmeta$cancertype <- rep(ct, nrow(tempmeta))
  # add to combined meta dataframe
  if (is.null(metaTCGA)){
    metaTCGA <- tempmeta
  } else{
    metaTCGA <- rbind(metaTCGA, tempmeta)
  }
  # extract cancer speicfic expression data
  tempexp <- TCGA_data[[ct]][["exp"]]
  tempexp$cancertype <- rep(ct, nrow(tempexp))
  # add to combined meta dataframe
  if (is.null(expTCGA)){
    expTCGA <- tempexp
  } else{
    expTCGA <- rbind(expTCGA, tempexp)
  }
}

#############################################################################
# save TCGA R data structures for future use
#############################################################################


#############################################################################
# Load CoGAPS results object and preprocessing TCGA data
#############################################################################
  ## load CoGAPS result object
  load("myeloid_run7_result3.RData")
  # restrict TCGA expression matrix to rows for which CoGAPS estimated feature
  # loadings
  #mtch <- match(toupper(rownames(gapsResult@featureLoadings)),
                #toupper(rownames(expTCGA)))
  #expTCGA <- expTCGA[mtch[!is.na(mtch)],]
  #gapsResult@featureLoadings <- gapsResult@featureLoadings[which(!is.na(mtch)),]
  symbols <- rownames(expTCGA)
  ###########################################################################
  # Run transfer learning
  ###########################################################################
  gaps2TCGA <- projectR(data = expTCGA, 
                      loadings = gapsResult,
                      full = TRUE, 
                      dataNames = symbols)
  ## get projected matrix out of new TL object
  TL.proj <- gaps2TCGA[["projection"]]
  ## set any NAs to zero
  TL.proj[which(is.na(TL.proj))] <- 0
  ## confirm NAs were removed and no other atypical value types (Inf, Nan, etc)
  stopifnot(!is.na(TL.proj))
  stopifnot(!is.infinite(TL.proj))
  stopifnot(!is.nan(TL.proj))

  ## save TL results REVISIT
  #saveRDS(TL.proj, file = "human_ici_scrna_projection.rds")
  # read in if starting analysis from here and don't want to rerun steps above
  #TL.proj <- readRDS("human_ici_scrna_projection.rds")

  
  # NOTE TO SELF: binary color coding to visualize enrichment
  
  ############################################################################
  # Graph results
  #############################################################################
  ct <- "TCGA-SKCM"
  View(metaTCGA)
  # Plot Complex Heatmaps
  # set cols - use viridis palette
  #mycols <- viridis(length(levels(simple.row_anno)), option = "D")
  mycols <- rainbow(6)
  # open connection to plotting device
  tiff(file = paste("TL.all_CoGAPs_patterns.ComplexHeatmap.TCGA.",ct,".tiff"),
       width = 6000, height = 5000, units = "px", res = 800)
  row_ha = rowAnnotation("Survival (months)" =
                           anno_lines(metaTCGA$days_to_death,
                                      axis_param= list(direction = "reverse",
                                                       labels_rot = 90),
                                      add_points = F))
                           # anno_lines(metaTCGA$paper_Survival..months.,
                           #            axis_param= list(direction = "reverse",
                           #                             labels_rot = 90),
                           #            add_points = F),
                         # "ImmuneScore" =
                         #   anno_lines(metaTCGA$paper_ESTIMATE.immune.score,
                         #              axis_param= list(direction = "reverse",
                         #                               labels_rot = 90),
                         #              add_points = F),
                         #   annotation_name_rot = 90)
  #scale data
  TL.proj.rownames <- rownames(TL.proj)
  TL.proj.scl <- apply(TL.proj, 2, scale)
  rownames(TL.proj.scl) <- TL.proj.rownames # reassign col as row names
  Heatmap(t(TL.proj.scl), col=inferno(100), name = "mat",
          clustering_distance_rows = "kendall",
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
Heatmap(t(TL.proj.pre), col=inferno(100), name = "mat",
        clustering_distance_rows = "kendall",
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
Heatmap(t(TL.proj.post), col=inferno(100), name = "mat",
        clustering_distance_rows = "kendall",
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

# transpose TL projection data and convert to df
TL.proj.df = as.data.frame(t(TL.proj))
TL.proj.df$group <- rownames(TL.proj.df)
TL.proj.df$status <- sfmeta$status
#mdat$cluster <- pdat$cluster_id
#mdat$response <- pdat$response
#mdat$combined <- pdat$combined
#mdat$group <- paste(pdat$combined, pdat$cluster_id, sep = "_")
TL.proj.df.m <- reshape2::melt(TL.proj.df, varnames = c("group", "status"),
             na.rm = FALSE, as.is = FALSE, value.name = "value")


colnames(TL.proj.df.m) <- c("group", "status", "pattern", "weight")
TL.proj.df.m$group <- gsub(".", " ", TL.proj.df.m$group, fixed=T)
TL.proj.df.m$group <- gsub("[[:digit:]]", "", TL.proj.df.m$group)
TL.proj.df.m$group <- gsub("CTLA", "CTLA4", TL.proj.df.m$group)
TL.proj.df.m$group <- gsub("ponder ", "ponder", TL.proj.df.m$group)
TL.proj.df.m$group <- as.factor(TL.proj.df.m$group)
# TL.proj.df.m$group <- factor(TL.proj.df.m$group, levels = c("Pre anti-CTLA4 Responder",
#                                             "Pre anti-CTLA4 Non-responder",
#                                             "Pre anti-PD1 Responder",
#                                             "Pre anti-PD1 Non-responder",
#                                             "Pre anti-CTLA4+PD1 Responder",
#                                             "Pre anti-CTLA4+PD1 Non-responder",
#                                             "Post anti-PD1 Responder",
#                                             "Post anti-PD1 Non-responder",
#                                             "Post anti-CTLA4+PD1 Responder",
#                                             "Post anti-CTLA4+PD1 Non-responder"))


## add partition data
# ask emily
# View(tl_cds)
# mdat$partition <- partitions(tl_cds)
# mdat$cluster_id <- pData(tl_cds)$cluster_id
## loop to plot all patterns by group 
# num_pats <- nrow(TL.proj)
# plot_list = list()
# 
# for(i in 1:num_pats){
#   pdat$test <- gapsResult@sampleFactors[,i]
#   plotname <- paste("Pattern", i, "of", num_pats, "in all cells", sep = " ")
#   p <- ggplot(pdat, aes(x=treatment, y=test, fill = treatment)) +
#     geom_violin(alpha = 0.8) + theme_bw() +
#     theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
#           axis.ticks.x=element_blank()) + 
#     scale_fill_viridis_d() + ggtitle(plotname) +
#     geom_point(aes(fill = treatment), size = 1.5, shape = 21, position = position_jitterdodge()) + 
#     ylab("pattern weight") + facet_wrap(~assigned_cell_type)
#   plot_list[[i]] = p
# }
# ## save plots together
# pdf("21og_pats_violin_assigned_celltype.pdf")
# for (i in 1:num_pats) {
#   print(plot_list[[i]])
# }
# dev.off()

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

## FIGURE 3B
pmdat <- mdat[mdat$pattern == num_pats[7],]
pmdat <- pmdat[pmdat$group == "Pre anti-CTLA4 Responder" | pmdat$group == "Pre anti-CTLA4 Non-responder",]
tiff(file = "tl_ctla4only_og.tiff", width = 6000, height = 5000, units = "px", res = 800)

ggplot(pmdat, aes(x=group, y=weight, fill = group)) +
  geom_violin(alpha = 0.8) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) + 
  scale_fill_viridis_d() + 
  geom_point(alpha = 0.8, aes(fill = group), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") 
dev.off()

t.test(pmdat[pmdat$group == "Pre anti-CTLA4 Responder",]$weight, pmdat[pmdat$group == "Pre anti-CTLA4 Non-responder",]$weight)


## FIGURE 3A
## subset for ctla4 monotherapy groups
pData(tl_cds)$test <- AP22[7,]
ctla4 <- tl_cds[,pData(tl_cds)$combined == "Pre anti-CTLA4 Responder" | pData(tl_cds)$combined == "Pre anti-CTLA4 Non-responder"]

ctla4 <- preprocess_cds(ctla4, num_dim = 5)
## remove batch effect of 
#ctla4 <- align_cds(ctla4, alignment_group = "patient")
plot_pc_variance_explained(ctla4)

## Reduce the dimensions using UMAP
ctla4 <- reduce_dimension(ctla4, umap.n_neighbors = 20L, max_components = 2)
## Cluster the cells
ctla4 <- cluster_cells(ctla4)
plot_cells(ctla4, color_cells_by = "combined")
plot_cells(ctla4, color_cells_by = "patient")

## plot UMAP by projected pattern weight
plot_cells(ctla4, color_cells_by = "test")

pdat <- data.frame(pData(ctla4))
pdat$umap1 <- ctla4@reducedDims@listData[["UMAP"]][,1]
pdat$umap2 <- ctla4@reducedDims@listData[["UMAP"]][,2]

cairo_ps("tl_ctla4_mono_umap_flat_p7.eps", width = 6, height = 6)   
tiff(file = "tl_ctla4_mono_umap_flat_p7.tiff", width = 6000, height = 5000, units = "px", res = 800)

ggplot(pdat, aes(x = umap1, y = umap2, color = test)) +
  geom_point(size=1.5) + 
  scale_color_viridis(option="magma", 
                      guide = guide_colorbar(barwidth = 1, barheight = 28, title="pattern weight",
                                             frame.colour = "black")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", 
                                        size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()


tiff(file = "tl_ctla4_mono_umap_flat_responder2.tiff", width = 6000, height = 5000, units = "px", res = 800)

#cairo_ps("tl_ctla4_mono_umap_flat_responder.eps", width = 6, height = 6)        
ggplot(pdat, aes(x = umap1, y = umap2, color = response)) +
  geom_point(size=1.5, aes(color=response)) + 
  scale_colour_viridis_d(direction = -1) + theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

###############################################################################
## 3D UMAP by response
pdat <- data.frame(pData(ctla4))
col_palette <- c("black", "purple")
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = ctla4@reducedDims@listData[["UMAP"]][,1], y = ctla4@reducedDims@listData[["UMAP"]][,2], z = ctla4@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$response)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$response)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("TL_3D_ctla4_mono_responders.png")

## 3D UMAP by pattern 7
par3d(windowRect = c(100, 100, 500,500))
plot3d(x = ctla4@reducedDims@listData[["UMAP"]][,1], y = ctla4@reducedDims@listData[["UMAP"]][,2], z = ctla4@reducedDims@listData[["UMAP"]][,3],
       col = myColorRamp(col = magma(n=10), values = AP22[7,]),
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
#par3d(pp)
rgl.snapshot("TL_3D_ctla4_mono_p7.png")

colkey (side = 3, add = TRUE, clim = myColorRamp(col = magma(n=10), values = AP22[7,]))


gapsnames <- rownames(gapsResult@featureLoadings)
subcounts <- counts2df[rownames(counts2df) %in% gapsnames,]
subcounts <- subcounts[,c(1:16290)]
## make cds to plot projected patterns 
expression_matrix <- as.matrix(subcounts)
gene_annotation <- data.frame(gene_short_name = rownames(subcounts))
rownames(gene_annotation) <- rownames(subcounts)
cell_metadata <- (meta2[meta2$title %in% colnames(subcounts),])
rownames(cell_metadata) <- colnames(expression_matrix)

tl_cds <- new_cell_data_set(expression_matrix,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)

save(tl_cds, file = "tl_cds_sadefeldman.rda")

load("tl_cds_sadefeldman.rda")

pData(tl_cds)$patient <- factor(pData(tl_cds)$patient)
tl_cds <- preprocess_cds(tl_cds, num_dim = 18)
## remove batch effect of 
tl_cds <- align_cds(tl_cds, alignment_group = "patient")
plot_pc_variance_explained(tl_cds)

## Step 2: Reduce the dimensions using UMAP
tl_cds <- reduce_dimension(tl_cds, umap.n_neighbors = 50L, max_components = 3)
tl_cds <- reduce_dimension(tl_cds, umap.n_neighbors = 200L, max_components = 2)
## Step 3: Cluster the cells
tl_cds <- cluster_cells(tl_cds)
plot_cells(tl_cds, color_cells_by = "combined")

marker_genes <- c("CD3D", "CD4", "CD8A", "FOXP3", "NCR1")
plot_cells(tl_cds, genes = marker_genes)
pData(tl_cds)$test <- AP22[7,]
plot_cells(tl_cds, color_cells_by = "test")

plot_cells(tl_cds, color_cells_by = "patient")
pData(tl_cds)$cluster_id <- factor(pData(tl_cds)$cluster_id)
plot_cells(tl_cds, color_cells_by = "cluster_id")
plot_cells(tl_cds, color_cells_by = "partition")

pData(tl_cds)$Cluster <- tl_cds@clusters@listData[["UMAP"]][["clusters"]]
pdat <- data.frame(pData(tl_cds))
pdat <- droplevels(pdat)

col_palette <- rainbow(n=11)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$cluster_id)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$cluster_id)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("tl_3D_cluster.png")

col_palette <- viridis(n=10)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$combined)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
#par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$combined)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("tl_3D_treatments.png")

pData(tl_cds)$response <- factor(pData(tl_cds)$response)
col_palette <- c("black", "purple")
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$response)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$response)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("tl_3D_responders.png")

col_palette <- viridis(n=3)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$therapy)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$therapy)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("tl_3D_therapy.png")

col_palette <- inferno(n=2)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$status)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$status)), pch = 16, col = col_palette, cex=1, inset=c(0.02))
#rgl.postscript("tcells_3D_treatments.eps","eps")
rgl.snapshot("tl_3D_status.png")



col_palette <- rainbow(n=32)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$patient)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
rgl.snapshot("tl_3D_patient.png")
## add in annotations from supplement
setwd("/Users/emily/Downloads")
annot <- read.csv("cell_annot.csv", header = TRUE)
pdat$cluster_id <- annot$Cluster.number[1:16290]

pdat$response <- droplevels(pdat$response)
pdat$therapy <- droplevels(pdat$therapy)

## plot flat umaps with scale legend
plot_cells(tl_cds)
pData(tl_cds)$test <- AP22[14,]
plot_cells(tl_cds, color_cells_by = "response")
plot_cells(tl_cds, color_cells_by = "therapy")
plot_cells(tl_cds, color_cells_by = "patient")
plot_cells(tl_cds, color_cells_by = "partition")

pdat <- data.frame(pData(tl_cds))
## add umap coordiantes to pdat
pdat$umap1 <- tl_cds@reducedDims@listData[["UMAP"]][,1]
pdat$umap2 <- tl_cds@reducedDims@listData[["UMAP"]][,2]
pdat$umap3 <- tl_cds@reducedDims@listData[["UMAP"]][,3]

pdat$test <- AP22[7,]
pdat$fill <- "x"

cairo_ps("umap_flat_tl_p7.eps", width = 6, height = 6)        
ggplot(pdat, aes(x = umap1, y = umap2, color = fill)) +
  geom_point(size=0.8) + 
  scale_color_manual(values = c("black")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())
dev.off()

ggplot(pdat, aes(x = umap1, y = umap2, color = combined)) +
  geom_point(size=0.5) + 
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

par3d(windowRect = c(100, 100, 500,500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col = myColorRamp(col = magma(n=10), values = AP22[7,]),
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
#par3d(pp)
#rgl.snapshot("tl_p5_21ogpats_3D.png")

colkey (side = 3, add = TRUE, clim = myColorRamp(col = magma(n=10), values = AP22[7,]))

legend3d("topright", legend = levels(as.factor(nkpdat$treatment)), pch = 16, col = treatment_palette, cex=1, inset=c(0.02))
rgl.snapshot("tl_p14_50pats_3D.png")

## add in sex info from supplement
pdat <- data.frame(pData(tl_cds))
#str_split_fixed(pdat$patient.ID..Pre.baseline..Post..on.treatment., "_", 2)
pdat$sex <- "F"
pdat[pdat$patient %in% c("P1","P2","P4","P5","P7","P8","P12","P13","P14","P15","P16","P18","P19",
                         "P22","P23","P24","P25","P26","P29","P30","P31","P35"),]$sex <- "M"

col_palette <- rainbow(n=5)
par3d(windowRect = c(100, 100, 500, 500))
plot3d(x = tl_cds@reducedDims@listData[["UMAP"]][,1], y = tl_cds@reducedDims@listData[["UMAP"]][,2], z = tl_cds@reducedDims@listData[["UMAP"]][,3],
       col=col_palette[as.factor(pdat$sex)],
       xlab="",ylab="",zlab="",size=2,alpha=.5,box=T,axes=T)
par3d(pp)
legend3d("topright", legend = levels(as.factor(pdat$sex)), pch = 16, col = treatment_palette, cex=1, inset=c(0.02))

rgl.snapshot("tl_3D_sex.png")

### TRAIN RANDOM FOREST CLASSIFIER 


library(caret)
## simplify pre/post and subset by therapy 
TL.proj.bin <- cbind(TL.proj)
colnames(TL.proj.bin) <- simple.row_anno
pd1 <- TL.proj.bin[,sfmeta$therapy == "anti-PD1"]
source("aucMat.R")
source("alluvialMat.R")
library(scales)
## check which patterns can predict which labels
## calculates AUC values
res <- aucMat(labels = colnames(TL.proj), weights = TL.proj)
res2 <- aucMat(labels = smeta$binary, weights = AP22)
View(t(res))
alluvialMat(gaps2mel, colnames(TL.proj), annotationName = "Treatment",
            annotationType = "Cell", plot = TRUE, minPropExplained = 0.9)

## 36 pd1/ctla4 responder
## 37 pd1 responders
## 21, 6 ctla4 responders
## 38 pd1 or combo responders
## 14 ctla4 or combo non responders

# We set a seed for the sake of reproducibility
set.seed(42)
class <- colnames(pd1)
## make colnames unique 
colnames(pd1) <- make.unique(colnames(pd1))
# test each pattern 
pat <- t(pd1)
rownames(pat) <- colnames(pd1)
# First, we'll pick off the training data: 
inTrain<-createDataPartition(y=class, p = 0.60, list=FALSE)
training<-(pat[inTrain,])
trainclass <- rownames(pat)[inTrain]
# Then what's left is testing data.
testing<-(pat[-inTrain,])
testclass <- rownames(pat)[-inTrain]

rf_model_1<-train(trainclass ~ ., data=training, method="rf")
confusionMatrix(predict(rf_model_1, training), reference=training$class, positive="1")



