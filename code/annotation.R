###################################################
### init
###################################################
library(monocle3)

###################################################
### code chunk: annotate cell types
###################################################
load(file = "monocle3_gubin_cds_init.rda")

# project the data onto the top principal components
cds <- preprocess_cds(cds,  method = 'PCA', 
                      norm_method = 'log', 
                      num_dim = 15, #28
                      verbose = T)
plot_pc_variance_explained(cds) 
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution = 5e-5, # 4e-3, 
                       reduction_method = 'UMAP')
plot_cells(cds, color_cells_by="cluster", label_cell_groups=TRUE)

## read in annotations 
## this file also contains the UMAP coordinates
anno <- read.csv(file = 'gubin_annotation.csv', header = TRUE, row.names = 1)
pdat <- data.frame(pData(cds))
pdat$assigned_cell_type <- ""
pdat[rownames(pdat) %in% rownames(anno),]$assigned_cell_type <- anno$assigned_cell_type
pData(cds)$assigned_cell_type <- pdat$assigned_cell_type
## remove unknowns
cds <- cds[,!(pData(cds)$assigned_cell_type == "")]

## to use the provided UMAP coordinates
cds@reducedDims@listData[["UMAP"]][,1] <- anno$umap1
cds@reducedDims@listData[["UMAP"]][,2] <- anno$umap2

plot_cells(cds, color_cells_by="assigned_cell_type", label_cell_groups=TRUE)

save(cds, file = "monocle3_gubin_cds_annotated.rda")

###################################################
### prepare matrix for CoGAPS
###################################################

p <- exprs(cds)
## remove genes with standard dev of 0
map <- fData(cds)
sum(apply(p,1,sd)==0,na.rm=T)
keepIndex = (apply(p,1,sd)!=0)
sum(keepIndex);sum(!keepIndex)
p = p[keepIndex,]
## log transform 
p <- log2(p+1)
## save sparse matrix
sparse.p <- Matrix(p, sparse = T)
## output for running CoGAPS 
writeMM(obj = sparse.p, file="gubin.mtx")

## this matrix was used as input for CoGAPS on AWS
