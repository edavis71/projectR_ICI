###################################################
# Code for plotting UMAPs and CoGAPS results for
# Gubin et al. 
# Author: Emily Davis-Marcisak
###################################################

###################################################
### init
###################################################
library(gplots)
library(ggplot2)
library(CoGAPS)
# load CoGAPs run result
load("data/myeloid_run7_result3.RData")
load("monocle3_gubin_cds_annotated.rda")
###################################################
### heatmap of all patterns
###################################################
pdat <- readRDS("gubin_annotation.rds")

## order factors for coloring
pdat$assigned_cell_type <- factor(pdat$assigned_cell_type, levels = c("NK", "CD8", "CD4", "Treg", "Mki67_hi", 
                                                                      "pDC", "DC", "Macrophages/Monocytes","Neutrophils"))
celltype_palette = c("#D53E4F", "#F46D43", "#FDAE61", "gold", "green4", "mediumpurple", "#41B6C4", "#225EA8", "#081D58")

gapsmat <- t(gapsResult@sampleFactors)
## heatmap of all patterns
cols <- celltype_palette

heatmap.2(gapsmat,col=inferno, 
          trace='none', # turns off trace lines inside the heat map
          #distfun=function(c) dist(c, method = "euclidean"),
          distfun=function(c) as.dist(1-cor(t(c))),
          cexCol=1.5,cexRow=0.3, scale = "column",
          hclustfun=function(x) hclust(x, method="average"),
          labCol = FALSE, # remove colnames
          colsep = 4936, dendrogram='none', key.title = NULL, keysize=1.2, key.par = list(cex=0.5),
          ColSideColors= cols[as.numeric(pdat$assigned_cell_type)])

###################################################
### plot UMAPs
###################################################

## cell type UMAP
ggplot(pdat, aes(x = umap1, y = umap2, color = assigned_cell_type)) +
  geom_point(size=0.5) + scale_color_manual(values=c("#D53E4F", "#F46D43", "#FDAE61", "gold", "green4", "mediumpurple", "#41B6C4", "#225EA8", "#081D58")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

## treatment UMAP
ggplot(pdat, aes(x = umap1, y = umap2, color = treatment)) +
  geom_point(size=0.5) + scale_color_manual(values=c("#000004FF", "#D95F02", "#21908CFF", "#7570B3"), 
                                            guide=guide_legend(override.aes = list(size = 6, shape = 15), title = NULL)) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(), 
        legend.key = element_rect(colour = 'black'))

## pattern 12
pdat$p12 <- gapsResult@sampleFactors[,12]    
ggplot(pdat, aes(x = umap1, y = umap2, color = p12)) +
  geom_point(size=1) + 
  scale_color_viridis(option="magma", 
                      guide = guide_colorbar(barwidth = 1, barheight = 28, title="",
                                             frame.colour = "black")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

## pattern 13
pdat$p13 <- gapsResult@sampleFactors[,13]    
ggplot(pdat, aes(x = umap1, y = umap2, color = p13)) +
  geom_point(size=1) + 
  scale_color_viridis(option="magma", 
                      guide = guide_colorbar(barwidth = 1, barheight = 28, title="",
                                             frame.colour = "black")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())


## pattern 7 UMAP
pdat$p7 <- gapsResult@sampleFactors[,7]    
ggplot(pdat, aes(x = umap1, y = umap2, color = p7)) +
  geom_point(size=1) + 
  scale_color_viridis(option="magma", 
                      guide = guide_colorbar(barwidth = 1, barheight = 28, title="",
                                             frame.colour = "black")) +
  theme_bw() + labs(x="UMAP1", y = "UMAP2") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(colour = "black", size=0.5),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank())

###################################################
### boxplot of pattern 12 & 13 in macrophages 
### by treatment
###################################################

mac_pdat <- pdat[pdat$assigned_cell_type == "Macrophages/Monocytes",]

ggplot(mac_pdat, aes(x=treatment, y=p12, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c( "#000004FF", "#D95F02","#21908CFF", "#7570B3")) + 
  geom_point(aes(), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") + coord_cartesian(ylim = c(0, 1.2))

ggplot(mac_pdat, aes(x=treatment, y=p13, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c( "#000004FF", "#D95F02","#21908CFF", "#7570B3")) + 
  geom_point(aes(), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") + coord_cartesian(ylim = c(0, 1.2))

###################################################
### boxplots of pattern 7 weight by celltype 
###################################################

ggplot(pdat, aes(x=assigned_cell_type, y=p7, fill = assigned_cell_type)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values=celltype_palette) +
  geom_point(aes(fill = assigned_cell_type), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0, 1.2))

###################################################
### boxplot of pattern 7 in NK cells by treatment
###################################################
nkpdat <- pdat[pdat$assigned_cell_type == "NK",]

ggplot(nkpdat, aes(x=treatment, y=p7, fill = treatment)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c( "#000004FF", "#D95F02","#21908CFF", "#7570B3")) + 
  geom_point(aes(), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") + coord_cartesian(ylim = c(0, 1.2))

