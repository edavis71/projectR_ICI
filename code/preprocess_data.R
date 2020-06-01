###################################################
# Code for preprocessing Gubin et al. 2018 data 
# prior to analysis with CoGAPS
# Author: Emily Davis-Marcisak
###################################################

###################################################
### init
###################################################
library(monocle3)
library(Matrix)
library(GEOquery)
library(plyr)

###################################################
### download Gubin et al. data from GEO
###################################################

filePaths = getGEOSuppFiles("GSE119352")
filePaths
tarF <- list.files(path = "./GSE119352/", pattern = "*.tar", full.names = TRUE)
tarF
untar(tarF, exdir = "./GSE119352/")
## get all the zip files
gzipF <- list.files(path = "./GSE119352/", pattern = "*.gz", full.names = TRUE)
gzipF
## unzip all your files
ldply(.data = gzipF, .fun = gunzip)

###################################################
### load data and format
###################################################

## read in matrix count data matrixmarket format
control_matrix <- readMM(file = './GSE119352/GSM3371684_Control_matrix.mtx')
class(control_matrix)
control_matrix <- as.matrix(control_matrix)
## check matrix
control_matrix[1:5,1:5]
## check dims
dim(control_matrix)
## read in gene names 
control_genes <- read.table(file = './GSE119352/GSM3371684_Control_genes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(control_genes)
## set gene short names to rownames of matrix
rownames(control_matrix) <- control_genes[,2]
## check new rownames
control_matrix[1:5,1:5]
## read in sample names
control_barcodes <- read.table(file = './GSE119352/GSM3371684_Control_barcodes.tsv', 
                               sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(control_barcodes)
## set barcodes to colnames of matrix
colnames(control_matrix) <- control_barcodes[,1]
## check new colnames
control_matrix[1:5,1:5]
## add unique identifier to colnames
colnames(control_matrix) <- paste(colnames(control_matrix), "control", sep = "_")
## pData for treatment
control_pdat <- data.frame("samples" = colnames(control_matrix), "treatment" = "control")
## need to do this for all treatments and merge matrices
## read in matrix count data
aPD1_matrix <- readMM(file = './GSE119352/GSM3371685_aPD1_matrix.mtx')
aPD1_matrix <- as.matrix(aPD1_matrix)
## check matrix
aPD1_matrix[1:5,1:5]
## check dims
dim(aPD1_matrix)
## read in gene names 
aPD1_genes <- read.table(file = './GSE119352/GSM3371685_aPD1_genes.tsv', 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(aPD1_genes)
## set gene short names to rownames of matrix
rownames(aPD1_matrix) <- aPD1_genes[,2]
## check new rownames
aPD1_matrix[1:5,1:5]
## read in sample names
aPD1_barcodes <- read.table(file = './GSE119352/GSM3371685_aPD1_barcodes.tsv', 
                            sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(aPD1_barcodes)
## set barcodes to colnames of matrix
colnames(aPD1_matrix) <- aPD1_barcodes[,1]
## check new colnames
aPD1_matrix[1:5,1:5]
## add unique identifier to colnames
colnames(aPD1_matrix) <- paste(colnames(aPD1_matrix), "aPD1", sep = "_")
## pData for treatment
aPD1_pdat <- data.frame("samples" = colnames(aPD1_matrix), "treatment" = "aPD1")
## read in matrix count data
aCTLA4_matrix <- readMM(file = './GSE119352/GSM3371686_aCTLA4_matrix.mtx')
aCTLA4_matrix <- as.matrix(aCTLA4_matrix)
## check matrix
aCTLA4_matrix[1:5,1:5]
## check dims
dim(aCTLA4_matrix)
## read in gene names 
aCTLA4_genes <- read.table(file = './GSE119352/GSM3371686_aCTLA4_genes.tsv', 
                           sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(aCTLA4_genes)
## set gene short names to rownames of matrix
rownames(aCTLA4_matrix) <- aCTLA4_genes[,2]
## check new rownames
aCTLA4_matrix[1:5,1:5]
## read in sample names
aCTLA4_barcodes <- read.table(file = './GSE119352/GSM3371686_aCTLA4_barcodes.tsv', 
                              sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(aCTLA4_barcodes)
## set barcodes to colnames of matrix
colnames(aCTLA4_matrix) <- aCTLA4_barcodes[,1]
## check new colnames
aCTLA4_matrix[1:5,1:5]
## add unique identifier to colnames
colnames(aCTLA4_matrix) <- paste(colnames(aCTLA4_matrix), "aCTLA4", sep = "_")
## pData for treatment
aCTLA4_pdat <- data.frame("samples" = colnames(aCTLA4_matrix), "treatment" = "aCTLA4")
## read in matrix count data
combo_matrix <- readMM(file = './GSE119352/GSM3371687_aPD1-aCTLA4_matrix.mtx')
combo_matrix <- as.matrix(combo_matrix)
## check matrix
combo_matrix[1:5,1:5]
## check dims
dim(combo_matrix)
## read in gene names 
combo_genes <- read.table(file = './GSE119352/GSM3371687_aPD1-aCTLA4_genes.tsv', 
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(combo_genes)
## set gene short names to rownames of matrix
rownames(combo_matrix) <- combo_genes[,2]
## check new rownames
combo_matrix[1:5,1:5]
## read in sample names
combo_barcodes <- read.table(file = './GSE119352/GSM3371687_aPD1-aCTLA4_barcodes.tsv', 
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)
## check dims
dim(combo_barcodes)
## set barcodes to colnames of matrix
colnames(combo_matrix) <- combo_barcodes[,1]
## check new colnames
combo_matrix[1:5,1:5]
## add unique identifier to colnames
colnames(combo_matrix) <- paste(colnames(combo_matrix), "combo", sep = "_")
## pData for treatment
combo_pdat <- data.frame("samples" = colnames(combo_matrix), "treatment" = "aPD1-aCTLA4")

###################################################
### join data to create cds
###################################################

## join matrices
joined <- cbind(control_matrix,aPD1_matrix,aCTLA4_matrix,combo_matrix)
## check that gene length stayed the same, samples merged
dim(joined)
## make pData dictionary
pdat <- rbind(control_pdat, aPD1_pdat, aCTLA4_pdat, combo_pdat)
## set samples to rownames
rownames(pdat) <- pdat$samples
## generate fData
fdat <- toupper(as.matrix(control_genes))
## set gene short names as rownames - note hack no longer works in R 3.5.1 and will force unique rownames
rownames(fdat) <- fdat[,2]
fdat <- data.frame(fdat)
common_colnames <- c("ensembl_id", "gene_short_name")
colnames(fdat) <- common_colnames
## set matrix gene names
rownames(joined) <- rownames(fdat)
## create cds of merged data 
cds <- new_cell_data_set(joined,
                  cell_metadata = pdat,
                  gene_metadata = fdat)

## check summary data 
table(pData(cds)$treatment)

save(cds, file = "monocle3_gubin_cds_init.rda")
