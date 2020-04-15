###################################################
### code chunk: init
###################################################
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
setwd('/Users/emily/OneDrive - Johns Hopkins University/projects/NK')

###################################################
### code chunk: load data from GEO - you only need to do this once
###################################################

filePaths = getGEOSuppFiles("GSE120575")
filePaths
## get all the zip files
gzipF <- list.files(path = "./GSE120575/", pattern = "*.gz", full.names = TRUE)
gzipF
## unzip all your files
ldply(.data = gzipF, .fun = gunzip)

###################################################
### code chunk: preprocessing human single cell dataset (sade-feldman)
###################################################

## load TPM data
sfcounts <- data.table::fread(file = "./GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", fill = TRUE)
dim(sfcounts)
## 55739 16293
## check out format
sfcounts[1:5,1:5]

## save barcodes to check metadata
barcodes <- sfcounts[c(1),]
## remove extra metadata rows
sfcounts <- sfcounts[-c(1:2),]
## check out format again, looking at end
sfcounts[1:5,16289:16293]
## remove extra column at end 
tmp <- sfcounts[,-c(16293)]
## convert to dataframe cause data table is weird
df <- data.frame(tmp)
## set gene names to rownames
rownames(df) <- df[,1]
## now remove extra name column
df <- df[,-c(1)]
## check 
df[1:5,1:5]
dim(df)

## preprocess meta data
sfmeta <- read.table(file = "./GSE120575/GSE120575_patient_ID_single_cells.txt", sep = '\t', header = FALSE)
dim(sfmeta)
## check end
sfmeta[16306:16343,]
## cut out extra  at end
tmp2 <- sfmeta[c(1:16306),]
## cut out extra rows at end
tmp2 <- tmp2[,-c(12:35)]
## cut first 13 rows
tmp2 <- tmp2[-c(1:15),]
## check
tmp2[1:5,1:5]

colnames(tmp2) <- c("Sample name", "title", "source name", "organism", "patient ID (Pre=baseline; Post= on treatment)", 
                     "response", "therapy")
head(tmp2)
## remove remaining empty cols
tmp2 <- tmp2[,-c(8:14)]

## set df colnames to barcode names
colnames(df) <- tmp2$title

## split and combine meta information for grouping
split <- strsplit(as.character(tmp2[,5]),split='_', fixed=TRUE)
p1 <- sapply(split, "[", 1)
p2 <- sapply(split, "[", 2)
tmp2$patient <- p2
tmp2$status <- p1
tmp2$combined <- paste(tmp2$status, tmp2$therapy, tmp2$response, sep = " ")

readmat <- as.matrix(df)
## convert gene names to the same format as the CoGAPS object
symbols <- toupper(rownames(df))

## checking summary table of treatments
table(tmp2$combined)

###################################################
### code chunk: transfer learning
###################################################
## load CoGAPS result object
load("data/myeloid_run7_result3.RData")

gaps2mel2 <- projectR(data = readmat, 
                      loadings = gapsResult,
                      full = TRUE, 
                      dataNames = symbols)

## get projected matrix out of new TL object
AP22 <- gaps2mel2[["projection"]]

## add back in sample or treatment names
colnames(AP22) <- tmp2$combined
tmp2 <- droplevels(tmp2)

## remove any NAs
AP22[which(is.na(AP22))] <- 0
## confirm NAs were removed
#which(is.na(AP22))
#which(is.infinite(AP22))
which(is.nan(AP22))

## save TL result so you don't have to do this again
#saveRDS(AP22, file = "human_ici_scrna_projection.rds")
AP22 <- readRDS("meta_human_ici_projection_ogrun.rds")
#saveRDS(AP22, file = "meta_human_ici_projection_ogrun.rds")
APog <- readRDS("meta_human_ici_projection_ogrun.rds")

## testing some heatmaps
graphics.off()
par(mar=c(7,4,4,2)+0.5) 

cols=brewer.pal(3,"Spectral")

tiff(file = "figures/heatmap_TL_human_combo", width = 6000, height = 5000, units = "px", res = 800)

row_ha = rowAnnotation(foo2 = as.factor(colnames(AP22)))

Heatmap(t(AP22), col=inferno(100), name = "mat", 
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 2,
        #show_column_names = FALSE
)


dev.off()

m = t(as.matrix(AP22))
mdat <-as.data.frame(m)
mdat$group <- ""
mdat$group <- rownames(m)
mdat$status <- pdat$status
#mdat$cluster <- pdat$cluster_id
#mdat$response <- pdat$response
#mdat$combined <- pdat$combined
#mdat$group <- paste(pdat$combined, pdat$cluster_id, sep = "_")
mdat <- melt(m, varnames = names(dimnames(m)),
             na.rm = FALSE, as.is = FALSE, value.name = "value")

colnames(mdat) <- c("group", "pattern", "weight")
mdat$group <- factor(mdat$group, levels = c("Pre anti-CTLA4 Responder",
                                            "Pre anti-CTLA4 Non-responder",
                                            "Pre anti-PD1 Responder",
                                            "Pre anti-PD1 Non-responder",
                                            "Pre anti-CTLA4+PD1 Responder",
                                            "Pre anti-CTLA4+PD1 Non-responder", 
                                            "Post anti-PD1 Responder",
                                            "Post anti-PD1 Non-responder",
                                            "Post anti-CTLA4+PD1 Responder",
                                            "Post anti-CTLA4+PD1 Non-responder"))


## add partition data
mdat$partition <- partitions(tl_cds)
mdat$cluster_id <- pData(tl_cds)$cluster_id
## loop to plot all patterns by group 
num_pats <- rownames(AP22)

plot_list = list()
for(i in 1:num_pats){
  pdat$test <- gapsResult@sampleFactors[,i]
  plotname <- paste("Pattern", i, "of", num_pats, "in all cells", sep = " ")
  p <- ggplot(pdat, aes(x=treatment, y=test, fill = treatment)) +
    geom_violin(alpha = 0.8) + theme_bw() +
    theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
          axis.ticks.x=element_blank()) + 
    scale_fill_viridis_d() + ggtitle(plotname) +
    geom_point(aes(fill = treatment), size = 1.5, shape = 21, position = position_jitterdodge()) + 
    ylab("pattern weight") + facet_wrap(~assigned_cell_type)
  plot_list[[i]] = p
}
## save plots together
pdf("21og_pats_violin_assigned_celltype.pdf")
for (i in 1:num_pats) {
  print(plot_list[[i]])
}
dev.off()

## violin plot of projected pattern weights
num_pats <- rownames(AP22)
plot_list = list()
for(i in 1:length(num_pats)){
  pmdat <- mdat[mdat$pattern == num_pats[i],]
  plotname <- paste("Pattern", i, "of", num_pats[i], sep = " ")
  p <- ggplot(pmdat, aes(x=group, y=weight, fill = group)) +
    geom_violin(alpha = 0.8) + theme_bw() +
    theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
          axis.ticks.x=element_blank()) + 
    scale_fill_viridis_d() + ggtitle(plotname) +
    geom_point(alpha = 0.8, aes(fill = group), size = 1.5, shape = 21, position = position_jitterdodge()) + 
    ylab("pattern weight") 
  plot_list[[i]] = p
}
## save plots together
pdf("21pats_TL_sc_og_violinplots_trim.pdf", width = 7, height = 6)
for (i in 1:length(num_pats)) {
  print(plot_list[[i]])
}
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
## removing pre/post 
meta2$binary <- paste(meta2$therapy, meta2$response, sep = " ")
## subset by therapy 
smeta <- meta2[1:16290,]
colnames(AP22) <- smeta$binary
pd1 <- AP22[,smeta$therapy == "anti-PD1"]

## check which patterns can predict which labels
## calculates AUC values
res <- aucMat(labels = smeta$combined, weights = AP22)
res2 <- aucMat(labels = smeta$binary, weights = AP22)
View(res)
alluvialMat(gaps2mel2, smeta$binary, annotationName = "Treatment",
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



