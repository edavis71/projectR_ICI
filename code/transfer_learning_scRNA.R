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
pdat_sf <- sfmeta[c(1:16306),]
## cut out extra cols at end
pdat_sf <- pdat_sf[,-c(12:35)]
## cut first 13 rows
pdat_sf <- pdat_sf[-c(1:15),]
## check
pdat_sf[1:5,1:5]

colnames(pdat_sf) <- c("Sample name", "title", "source name", "organism", "patient ID (Pre=baseline; Post= on treatment)", 
                    "response", "therapy")
head(pdat_sf)
## remove remaining empty cols
pdat_sf <- pdat_sf[,-c(8:14)]
## set df colnames to barcode names
colnames(df) <- pdat_sf$title

## split and combine meta information for grouping
split <- strsplit(as.character(pdat_sf[,5]),split='_', fixed=TRUE)
p1 <- sapply(split, "[", 1)
p2 <- sapply(split, "[", 2)
pdat_sf$patient <- p2
pdat_sf$status <- p1
pdat_sf$combined <- paste(pdat_sf$status, pdat_sf$therapy, pdat_sf$response, sep = " ")

## convert to numeric matrix
readmat <- data.matrix(df)
## convert gene names to the same format as the CoGAPS object
symbols <- toupper(rownames(df))

## summary table of treatments
table(pdat_sf$combined)

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
colnames(sf_proj) <- pdat_sf$combined
pdat_sf <- droplevels(pdat_sf)

saveRDS(sf_proj, file = "sf_projection.rds")

###################################################
### plotting transfer learning: sade-feldman 
### all cells
###################################################

## add pattern 7 weight to meta data
pdat_sf$p7 <- sf_proj[7,]

sf_pre <- pdat_sf[pdat_sf$status == "Pre",]
## pattern 7 weight across all cells
ggplot(sf_pre, aes(x=combined, y=p7, color= combined, fill = combined)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("Pre anti-CTLA4 Responder" = "black", 
                                 "Pre anti-CTLA4 Non-responder" = "#21908CFF",
                                 "Pre anti-PD1 Non-responder" = "#D95F02", 
                                 "Pre anti-PD1 Responder" = "black",
                                 "Pre anti-CTLA4+PD1 Responder" =  "black",
                                 "Pre anti-CTLA4+PD1 Non-responder" = "#7570B3")) +
  scale_fill_manual(values = c("Pre anti-CTLA4 Responder" = "#21908CFF", 
                               "Pre anti-CTLA4 Non-responder" = "white",
                               "Pre anti-PD1 Non-responder" = "white", 
                               "Pre anti-PD1 Responder" = "#D95F02",
                               "Pre anti-CTLA4+PD1 Responder" =  "#7570B3",
                               "Pre anti-CTLA4+PD1 Non-responder" = "white")) + 
  geom_point(aes(color =combined), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  ylab("pattern weight") + xlab("treatment") + 
  facet_wrap(~therapy, scales = "free_x") + coord_cartesian(ylim = c(0, 2.5))

sf_post <- pdat_sf[pdat_sf$status == "Post",]
ggplot(sf_post, aes(x=combined, y=p7, color= combined, fill = combined)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("Post anti-PD1 Non-responder" = "#D95F02", 
                                 "Post anti-PD1 Responder" = "black",
                                 "Post anti-CTLA4+PD1 Responder" =  "black",
                                 "Post anti-CTLA4+PD1 Non-responder" = "#7570B3")) +
  scale_fill_manual(values = c("Post anti-PD1 Non-responder" = "white", 
                               "Post anti-PD1 Responder" = "#D95F02",
                               "Post anti-CTLA4+PD1 Responder" =  "#7570B3",
                               "Post anti-CTLA4+PD1 Non-responder" = "white")) + 
  geom_point(aes(color =combined), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") +
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  facet_wrap(~therapy, scales = "free_x") + coord_cartesian(ylim = c(0, 2.5))

###################################################
### plotting transfer learning: sade-feldman 
### NK cells
###################################################

t_mat <- as.matrix(t(readmat))
nk_gate <- c("NKG7", "FCGR3A", "NCR1", "CD3G", "CD3D", "CD4")
nk_mat <- t_mat[,colnames(t_mat) %in% nk_gate]
df <- as.data.frame(nk_mat)
df$nk <- ""
df[(df$NKG7 > 0 | df$FCGR3A > 0 | df$NCR1 > 0) &
   (df$CD3D == 0 & df$CD3G == 0 & df$CD4 == 0),]$nk <- "TRUE"

pdat_sf <- cbind(pdat_sf, df)
sf_nk <-  pdat_sf[pdat_sf$nk == "TRUE",]

sf_nk$combined <- factor(sf_nk$combined, levels = c("Pre anti-CTLA4 Non-responder",
                                                "Pre anti-CTLA4 Responder",
                                                "Pre anti-PD1 Non-responder",
                                                "Pre anti-PD1 Responder",
                                                "Pre anti-CTLA4+PD1 Non-responder",
                                                "Pre anti-CTLA4+PD1 Responder",
                                                "Post anti-PD1 Non-responder",
                                                "Post anti-PD1 Responder",
                                                "Post anti-CTLA4+PD1 Non-responder",
                                                "Post anti-CTLA4+PD1 Responder"))

sf_nk$therapy <- factor(sf_nk$therapy, levels = c("anti-CTLA4", "anti-PD1", "anti-CTLA4+PD1"))

## pre-treatment
sf_nk_pre <- sf_nk[sf_nk$status == "Pre",]
ggplot(sf_nk_pre, aes(x=combined, y=p7, color= combined, fill = combined)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("Pre anti-CTLA4 Responder" = "black", 
                                 "Pre anti-CTLA4 Non-responder" = "#21908CFF",
                                 "Pre anti-PD1 Non-responder" = "#D95F02", 
                                 "Pre anti-PD1 Responder" = "black",
                                 "Pre anti-CTLA4+PD1 Responder" =  "black",
                                 "Pre anti-CTLA4+PD1 Non-responder" = "#7570B3")) +
  scale_fill_manual(values = c("Pre anti-CTLA4 Responder" = "#21908CFF", 
                               "Pre anti-CTLA4 Non-responder" = "white",
                               "Pre anti-PD1 Non-responder" = "white", 
                               "Pre anti-PD1 Responder" = "#D95F02",
                               "Pre anti-CTLA4+PD1 Responder" =  "#7570B3",
                               "Pre anti-CTLA4+PD1 Non-responder" = "white")) + 
  geom_point(aes(color =combined), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  ylab("pattern weight") + xlab("treatment") + 
  facet_wrap(~therapy, scales = "free_x") + coord_cartesian(ylim = c(0, 2.5))

## on-treatment
sf_nk_post <- sf_nk[sf_nk$status == "Post",]
ggplot(sf_nk_post, aes(x=combined, y=p7, color= combined, fill = combined)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("Post anti-PD1 Non-responder" = "#D95F02", 
                                 "Post anti-PD1 Responder" = "black",
                                 "Post anti-CTLA4+PD1 Responder" =  "black",
                                 "Post anti-CTLA4+PD1 Non-responder" = "#7570B3")) +
  scale_fill_manual(values = c("Post anti-PD1 Non-responder" = "white", 
                               "Post anti-PD1 Responder" = "#D95F02",
                               "Post anti-CTLA4+PD1 Responder" =  "#7570B3",
                               "Post anti-CTLA4+PD1 Non-responder" = "white")) + 
  geom_point(aes(color =combined), size = 1.5, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") +
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  facet_wrap(~therapy, scales = "free_x") + coord_cartesian(ylim = c(0, 2.5))

###################################################
### ROC analysis: sade-feldman 
###################################################
tmp2$p7 <- sf_proj[7,]
## clta-4 model 
pdat$
ctla4 <- tmp2[which(tmp2$therapy == "anti-CTLA4"),]
ctla4 <- droplevels(ctla4)
## binarized response 
ctla4$binary <- factor(1*(ctla4$response == 'Responder'))
table(ctla4$binary)
table(ctla4$combined)
mylogit <- glm(binary ~ p7, data = ctla4, family = "binomial")
summary(mylogit)
prob=predict(mylogit,type=c("response"))
ctla4$prob=prob
pred <- prediction(prob, ctla4$binary)    
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc <- performance(pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc
## 0.7482411

plot(perf, col="#21908CFF", lwd = 3)
abline(0, 1) #add a 45 degree line

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
### transfer learning: de Andrade
###################################################

combo_mat <- as.matrix(combo)
nk_symbols <- rownames(combo)
gaps2nk <- projectR(data = combo_mat, 
                    loadings = gapsResult,
                    full = TRUE, 
                    dataNames = nk_symbols)

## get projected matrix out of new TL object
nk_proj <- gaps2nk[["projection"]]
## add back in sample or treatment names
colnames(nk_proj) <- meta$patient

saveRDS(nk_proj, file = "da_projection.rds")

###################################################
### plotting transfer learning: de Andrade
###################################################

## fix labels
meta[meta$tissue == "tumor",]$tissue <- "Tumor"
meta[meta$tissue == "blood",]$tissue <- "Blood"
meta[meta$tissue == "Center",]$tissue <- "Tumor"
meta[meta$tissue == "Cortex",]$tissue <- "Tumor"
meta[meta$tissue == "Nodule",]$tissue <- "Tumor"

## add pattern 7 weight to meta data
meta$p7 <- nk_proj[7,]

## add in treatment annotations from supplemental file
meta$therapy <- ""
meta$resistance <- ""
meta$therapy = dplyr::recode(meta$patient,
                             "CY129" = "nivolumumab+pembrolizumab",
                             "CY155" = "pembrolizumab",
                             "CY158" = "ipilimumab",
                             "CY160" = "none",
                             "CY164" = "pembrolizumab+TVEC")
meta$resistance = dplyr::recode(meta$patient,
                                "CY129" = "primary",
                                "CY155" = "primary",
                                "CY158" = "acquired",
                                "CY160" = "primary", 
                                "CY164" = "acquired")

meta$patient <- factor(meta$patient, levels = c("CY160", "CY155", "CY164",
                                                "CY158", "CY129"))
tumor <- meta[meta$tissue == "Tumor",]
ggplot(tumor, aes(x=patient, y=p7, color= patient, fill = patient)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("CY160" = "#000004FF",
                                 "CY158" = "black",
                                 "CY155" = "#D95F02", 
                                 "CY164" = "black",
                                 "CY129" = "#7570B3")) +
  scale_fill_manual(values = c("CY129" = "white",
                               "CY155" = "white",
                               "CY158" = "#21908CFF",
                               "CY160" = "white", 
                               "CY164" = "#D95F02")) +  
  geom_point(aes(color =patient), size = 3, shape = 21, position = position_jitterdodge()) + 
  ylab("pattern weight") + xlab("treatment") + 
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  facet_wrap(~resistance, scales = "free_x") + coord_cartesian(ylim = c(-8.5, 18)) +
  geom_hline(yintercept=3.53, linetype="dashed", color = "black", size=0.5)

blood <- meta[meta$tissue == "Blood",]
ggplot(blood, aes(x=patient, y=p7, color= patient, fill = patient)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  scale_colour_manual(values = c("CY160" = "#000004FF",
                                 "CY158" = "black",
                                 "CY155" = "#D95F02", 
                                 "CY164" = "black",
                                 "CY129" = "#7570B3")) +
  scale_fill_manual(values = c("CY129" = "white",
                               "CY155" = "white",
                               "CY158" = "#21908CFF",
                               "CY160" = "white", 
                               "CY164" = "#D95F02")) +  
  geom_point(aes(color =patient), size = 3, shape = 21, position = position_jitterdodge()) + 
  theme(panel.spacing.x=unit(0.15, "lines"),panel.spacing.y=unit(1, "lines")) +
  ylab("pattern weight") + xlab("treatment") + 
  facet_wrap(~resistance, scales = "free_x") + coord_cartesian(ylim = c(-8.5, 18)) +
  geom_hline(yintercept=3.53, linetype="dashed", color = "black", size=0.5)


