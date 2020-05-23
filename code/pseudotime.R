###################################################
# Code for plotting pseudotime results for NK cells
# treated with aCTLA-4 
# Author: Emily Davis-Marcisak
###################################################

###################################################
### init
###################################################
library(monocle3)
library(gplots)
library(viridis)

load("data/myeloid_run7_result3.RData")
load("monocle3_myeloid_cds_updated.rda")
###################################################
### pseudotime
###################################################
pData(cds)$p7 <- gapsResult@sampleFactors[,7]  
nk <- cds[,pData(cds)$assigned_cell_type == "NK" & pData(cds)$treatment == "aCTLA4"]
nk <- preprocess_cds(nk, num_dim = 50)
plot_pc_variance_explained(nk)
nk <- reduce_dimension(nk, umap.n_neighbors = 23L, max_components = 2)
nk <- cluster_cells(nk)
plot_cells(nk, color_cells_by="p7", cell_size = 1.0)

## 50, 15/2
## trajectory analysis
nk <- learn_graph(nk)
nk <- order_cells(nk)

plot_cells(nk,
           color_cells_by = "p7",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size = 2.0)

plot_cells(nk,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.8)


pData(nk)$pseudo <- nk@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
cor.test(pData(nk)$pseudo, pData(nk)$p7)
###################################################
### differential expression
###################################################

pM <- patternMarkers(gapsResult)

pData(nk)$weight <- ""
pData(nk)[pData(nk)$p7 >= 0.35,]$weight <- "positive"
pData(nk)[pData(nk)$p7 < 0.35,]$weight <- "negative"

## calculate fold change
low_count <- sum(pData(nk)$weight == "negative" , na.rm=TRUE)
high_count <- sum(pData(nk)$weight == "positive", na.rm=TRUE)
fData(nk)$Total_mRNAs_per_gene_nk <- Matrix::rowSums(exprs(nk))
fData(nk)$Average_mRNAs_low <- Matrix::rowSums(exprs(nk[,pData(nk)$weight=="negative"]))/low_count
fData(nk)$Average_mRNAs_high <- Matrix::rowSums(exprs(nk[,pData(nk)$weight=="positive"]))/high_count
fData(nk)$Expression_change_nk <- ifelse(fData(nk)$Average_mRNAs_low < fData(nk)$Average_mRNAs_high, "Upregulated", ifelse(fData(nk)$Average_mRNAs_low > fData(nk)$Average_mRNAs_high, "Downregulated", "Unchanged"))
fData(nk)$Fold_change_nk <- fData(nk)$Average_mRNAs_high/fData(nk)$Average_mRNAs_low

## nkset for pattern markers
p7_genes = pM[["PatternMarkers"]][[7]]
p7_nkset = nk[rowData(nk)$gene_short_name %in% p7_genes,]
p7_fdat <- as.data.frame(fData(p7_nkset))

## identify genes that change as a function of pseudotime
test_res <- graph_test(nk, neighbor_graph="principal_graph", cores=4)
pr_deg <- subset(test_res, q_value < 0.01)

pr_deg_ids <- row.names(subset(test_res, q_value < 0.01))
## collect the trajectory-variable genes into modules: c(0,10^seq(-6,-1))
gene_module_df <- find_gene_modules(nk[pr_deg_ids,], resolution=1e-02)
cell_group_df <- tibble::tibble(cell=row.names(colData(nk)), 
                                cell_group=colData(nk)$bin)
agg_mat <- aggregate_gene_expression(nk, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

AFD_genes <- c("LY6A", "GZMF", "PRF1")
AFD_lineage_cds <- nk[rowData(nk)$gene_short_name %in% AFD_genes]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="p7",
                         min_expr=0.3)


## get list of top 100 most significant genes
mygenes <- pr_deg[order(pr_deg$q_value),]
## intersect with pattern markers
pm_sig <- mygenes[mygenes$gene_short_name %in% p7_fdat$gene_short_name,]
pm_sig <- pm_sig[pm_sig$Expression_change_nk == "Upregulated",]
pm_sig <- pm_sig[which(pm_sig$Fold_change_nk > 5.2),]

mat <- as.matrix(exprs(nk))
mat <- mat[rownames(mat) %in% pm_sig$gene_short_name,]
dim(mat)
colnames(mat) <- pData(nk)$p7

cols= plasma(n = 340)
#pal <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(log2(mat+1),col=inferno, 
          trace='none', # turns off trace lines inside the heat map
          distfun=function(c) dist(c, method = "euclidean"),
          #distfun=function(c) as.dist(1-cor(t(c))),
          cexCol=1.5,cexRow=0.3, scale = "row",
          hclustfun=function(x) hclust(x, method="average"),
          labCol = FALSE, # remove colnames
          ColSideColors= cols[as.numeric(as.factor(colnames(mat)))])
dev.off()

