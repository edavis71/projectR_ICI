###################################################
### init
###################################################
library(projectR)
library(CePa)
library(org.Hs.eg.db)
library(biomaRt)
library(ggplot2)
library(CoGAPS)
library(data.table)
library(stringr)
library(Matrix)
library(GEOquery)

# load CoGAPs result
load("data/gapsResult.RData")
###################################################
### load human scRNA-seq datasets GEO
###################################################

## sade-feldman et al. 2018
filePaths = getGEOSuppFiles("GSE120575")
filePaths
## de Andrade et al. 2019
filePaths2 = getGEOSuppFiles("GSE139249")
filePaths2

filePaths2 = getGEOSuppFiles("GSE139249")
filePaths2

###################################################
### preprocess sade-feldman data for TL
###################################################

## load TPM data
sfcounts <- data.table::fread(file = "GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", fill = TRUE)
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
## cut out extra rows at end
tmp2 <- sfmeta[c(1:16306),]
## cut out extra cols at end
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

## convert to numeric matrix
readmat <- data.matrix(df)
## convert gene names to the same format as the CoGAPS object
symbols <- toupper(rownames(df))

## summary table of treatments
table(tmp2$combined)

###################################################
### transfer learning: sade-feldman
###################################################

gaps2sf <- projectR(data = readmat, 
                     loadings = gapsResult,
                     full = TRUE, 
                     dataNames = symbols)

## get projected matrix out of new TL object
sf_proj <- gaps2sf[["projection"]]

## add back in sample or treatment names
colnames(sf_proj) <- tmp2$combined
tmp2 <- droplevels(tmp2)

saveRDS(sf_proj, file = "sf_projection.rds")
###################################################
### preprocess de Andrade data for TL
###################################################

setwd('./GSE139249')
dir.create('files')
untar("GSE139249_RAW.tar", exdir = "files")

## false returns a data frame
options(datatable.fread.datatable=FALSE)

setwd('./files')
files <- list.files()
split <- strsplit(files,split='_', fixed=TRUE)
gsm <- sapply(split, "[", 1)
gsm <- unique(gsm)

df_test <- list()
list_mtx <- list()
for (prefix in 1:length(gsm)) {
  gsm_name = gsm[prefix]
  print(gsm_name)
  
  expr_barcodes = paste(gsm_name, ".*", "barcodes.tsv.gz", sep = "")
  filename_barcodes <- files[grep(expr_barcodes, files)]
  df_barcodes <- data.table::fread(file = filename_barcodes, fill = TRUE, header = FALSE)
  sample <- str_extract(filename_barcodes, "CY[0-9]{3}")
  tissue <- rev(unlist(strsplit(filename_barcodes, "-|_")))[2]
  meta <- data.frame(barcodes = df_barcodes$V1, 
                     patient = sample, tissue = tissue)
  sample_name <- paste(df_barcodes$V1, rownames(meta), sep = "_")
  meta$sample_name <- sample_name
  df_test[[prefix]] <- meta
  
  expr_genes = paste(gsm_name, ".*", "genes.tsv.gz", sep = "")
  print(expr_genes)
  filename_genes <- files[grep(expr_genes, files)]
  df_genes <- data.table::fread(file = filename_genes, fill = TRUE, header = FALSE)
  print(dim(df_genes))
  expr_mtx = paste(gsm_name, ".*", "matrix.mtx.gz", sep = "")
  print(expr_mtx)
  filename_mtx <- files[grep(expr_mtx, files)]
  list_mtx[[prefix]] <- as.matrix(readMM(file = filename_mtx))
  rownames(list_mtx[[prefix]]) <- df_genes$V2
  colnames(list_mtx[[prefix]]) <- sample_name
}
## 33538 or 33694 genes, need to take only the ones that overlap
meta <- do.call(rbind, df_test)

combo <- list_mtx[[1]]
for (i in 2:length(list_mtx)) {
  combo <- merge(combo, list_mtx[[i]], by=0, all=F)
  rownames(combo) <- combo$Row.names
  combo <- combo[,-c(1)]
  print(dim(combo))
}

###################################################
### transfer learning: sade-feldman
###################################################

combo_mat <- as.matrix(combo)
nk.symbols <- rownames(combo)
gaps2nk <- projectR(data = combo_mat, 
                    loadings = gapsResult,
                    full = TRUE, 
                    dataNames = nk.symbols)

## get projected matrix out of new TL object
nk_proj <- gaps2nk[["projection"]]
## add back in sample or treatment names
colnames(nk_proj) <- meta$patient

saveRDS(nk_proj, file = "da_projection.rds")

