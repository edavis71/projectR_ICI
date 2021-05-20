#############################################################################
# Use the TCGAbiolinks R package to download all TCGA expression data
# Project CoGAPS patterns into TCGA data
# Perform survival analysis 
# Authors: Emily Davis-Marcisak
#############################################################################

library(SummarizedExperiment)
library(TCGAbiolinks)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(plyr)
library(tidyverse)
library(readxl)
library(CoGAPS)
library(projectR)
library(survival)

#############################################################################
# Using TCGAbiolinks to download data from GDC
#############################################################################

project_IDs <- TCGAbiolinks:::getGDCprojects()$project_id
#-- Get list of all TCGA named projects 
cancer_types <- project_IDs[grep("TCGA", TCGAbiolinks:::getGDCprojects()$project_id)]

tcga_list <- list() 
for (cancer in cancer_types){
  
  #-- Download legacy data for TCGA tumor types, Hg19
  query <- GDCquery(cancer, 
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  
  #-- Download harmonized data for TCGA tumor types, Hg38
  #query <- GDCquery(cancer, 
  #                  data.category = "Transcriptome Profiling", 
  #                  data.type = "Gene Expression Quantification",
  #                  workflow.type = "HTSeq - Counts") #workflow.type = "HTSeq - FPKM") 
  
  GDCdownload(query)
  #-- parse/prepare the data for each cancer type
  data <- GDCprepare(query)
  #--  add tumor type column to meta data
  SummarizedExperiment::colData(data)$cancer <- cancer
  #-- Subset by known gene locations
  geneRanges <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  data <- subsetByOverlaps(data, geneRanges)
  #-- store summarized experiment object in list
  tcga_list[[cancer]] <- data
}

save(tcga_list, file = "TCGA_legacy.rda")

#############################################################################
# Download survival data from the TCGA Pan-Cancer Clinical Data paper (Liu et. al, 2018)
# Join clinical features to SummarizedExperiment
#############################################################################

#-- Get more complete survival data (Liu et. al. 2018)
download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx",
              destfile = "Liu2018_survivalData.xlsx")

#-- Read survival data
nms <- names(read_excel("Liu2018_survivalData.xlsx", n_max = 0))
col_types <- c("skip", ifelse(grepl("residual_tumor|cause_of_death", nms), "text", "guess"))
liu_survData <- read_excel("Liu2018_survivalData.xlsx", sheet = 1, 
                           col_types = col_types, na = c("#N/A", "[Not Available]"))
liu_survData$cancer <- paste("TCGA-", liu_survData$type, sep ="")

for (cancer in cancer_types){
  #-- Get patient metadata for additional pre-processing
  colAnnotation <- colData(tcga_list[[cancer]])
  
  #-- Preprocess
  survData <- liu_survData %>%
    filter(cancer == cancer, 
           bcr_patient_barcode %in% colAnnotation$patient) %>%
    select(bcr_patient_barcode, gender, Age = age_at_initial_pathologic_diagnosis,
           ajcc_pathologic_tumor_stage, OS, OS.time, PFI, PFI.time) %>%
    mutate(OS.time.months = OS.time/30,
           PFI.time.months = PFI.time/30,
           Tumor_Stage = as.numeric(mapvalues(ajcc_pathologic_tumor_stage, 
                                              from = c("Stage I", "Stage II", "Stage III",
                                                       "Stage IIIA", "Stage IIIB", "Stage IIIC", 
                                                       "Stage IV", "Stage IVA", "Stage IVB",
                                                       "[Discrepancy]"),
                                              to = c(1, 2, 3, 3, 3, 3, 4, 4, 4, NA))),
          Stage = as.factor(mapvalues(Tumor_Stage,
                                       from = c(1,2,3,4), to = c("I", "II", "III", "IV")))) %>%
    as.data.frame()
    
    #-- join features to the data from GDC through TCGAbiolinks
    #-- Conform names
    rownames(survData) <- survData$bcr_patient_barcode
    survData <- survData[colAnnotation$patient,]
    rownames(survData) <- rownames(colAnnotation)
  
    #-- Add back to Summarized Experiment
    colData(tcga_list[[cancer]]) <- cbind(colData(tcga_list[[cancer]]), as(survData, "DataFrame"))
}
  
save(tcga_list, file = "TCGA_legacy_preprocessed.rda")

#############################################################################
# Transfer learning 21 mouse CoGAPS patterns to TCGA
#############################################################################

#-- Load CoGAPS result object for Gubin et al. 
load("data/myeloid_run7_result3.RData")

for (cancer in cancer_types){
  #-- Get expression data
  exp <- as.matrix(SummarizedExperiment::assay(tcga_list[[cancer]]))
  #-- Get gene symbols for projectR
  rowdata <- data.frame(rowData(tcga_list[[cancer]]))
  symbols <- rowdata$external_gene_name

  #-- Transfer learning
  gaps2tcga <- projectR(data = exp,
                        loadings = gapsResult,
                        full = TRUE,
                        dataNames = symbols)
  projection <- gaps2tcga[["projection"]]

  #-- Add back to SummarizedExperiment
  colData(tcga_list[[cancer]]) <- cbind(colData(tcga_list[[cancer]]), t(projection))
}

save(tcga_list, file = "TCGA_legacy_projection_gubin.rda")

#############################################################################
# Preprocess data for survival analysis
#############################################################################

meta_list <- list() 
for (cancer in cancer_types){
  df <- data.frame(colData(tcga_list[[cancer]]))
  #-- Store metadata for each single cell experiment 
  meta_list[[cancer]] <- df
}
#-- Subset dataframes to shared column names 
common_cols <- Reduce(intersect, lapply(meta_list, names))
meta_list <- lapply(meta_list, function(x) x[common_cols])
#-- Merge metadata for single cell experiment objects
meta_tcga <- do.call(rbind, meta_list)

#-- Combine definitions
meta_tcga$category <- meta_tcga$definition
meta_tcga[meta_tcga$definition %in% c("Primary solid Tumor", "Additional - New Primary"),]$category <- "Primary"
meta_tcga[meta_tcga$definition %in% c("Metastatic", "Additional Metastatic"),]$category <- "Metastatic"

save(meta_tcga,cancer_types,patterns, file = "TCGA_legacy_projection_gubin_df.rda")

#############################################################################
# Survival analysis 
# Cox Proportional Hazard Models for each TCGA cancer type and CoGAPS pattern 
#############################################################################

load("TCGA_legacy_projection_gubin_df.rda")
#-- Number of samples for each cancer type
table(meta_tcga$cancer)
table(meta_tcga$cancer, meta_tcga$category)
#-- Tabulate by outcome
xtabs(~cancer+OS, data=meta_tcga) %>% addmargins()

#-- Remove normal samples and non-solid tumors for analysis 
remove_list <- c("Solid Tissue Normal", "Primary Blood Derived Cancer - Peripheral Blood",
                 "Recurrent Solid Tumor", "Additional - New Primary", "Additional Metastatic")
no_norm <- meta_tcga[!meta_tcga$definition %in% remove_list,]
no_norm <- droplevels(no_norm)

table(no_norm$cancer, no_norm$category)

#-- Model cox survival for all cancer types across each pattern

#-- Univariate cox regression analysis 
#-- Apply univariate coxph function to multiple covariates 
univ_coxph <- function(covariates, data) {
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
  print(univ_formulas) ## CHECKING FORMULAS
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
  #-- Create dataframe of output 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value <- signif(x$wald["pvalue"], digits=2)
                           wald.test <- signif(x$wald["test"], digits=2)
                           beta <- signif(x$coef[1], digits=2); # coefficient beta
                           HR <- signif(x$coef[2], digits=2); # exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"], 2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res <- c(beta, HR, wald.test, p.value)
                           names(res) <- c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                         })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  as.data.frame(res)
  return(res)
}

covariates <- c("Pattern_1","Pattern_2","Pattern_3",
                "Pattern_4","Pattern_5","Pattern_6","Pattern_7",
                "Pattern_8","Pattern_9","Pattern_10","Pattern_11",
                "Pattern_12","Pattern_13","Pattern_14","Pattern_15",
                "Pattern_16","Pattern_17","Pattern_18","Pattern_19",
                "Pattern_20","Pattern_21","Age")

test <- univ_coxph(covariates, no_norm) 

#-- Multivariate cox regression analysis
#-- Age as covariates with each pattern
multiv_coxph <- function(covariates, data) {
  multiv_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(OS.time, OS)~Age+', x)))
  multiv_models <- lapply(multiv_formulas, function(x){coxph(x, data = data)})
  #return(multiv_models)
  #-- Create dataframe of output 
  multiv_results <- lapply(multiv_models,
                           function(x){ 
                             x <- summary(x)
                             pattern_coef <- x$coefficients[2,] # note that this will change depending on number of covariates included
                             res <- t(data.frame(pattern_coef))
                             rownames(res) <- rownames(x$coefficients)[2]
                             return(res)
                           })
  
  res <- do.call(rbind.data.frame, multiv_results)
  return(res)
}


covariates <- c("Pattern_1","Pattern_2","Pattern_3",
                "Pattern_4","Pattern_5","Pattern_6","Pattern_7",
                "Pattern_8","Pattern_9","Pattern_10","Pattern_11",
                "Pattern_12","Pattern_13","Pattern_14","Pattern_15",
                "Pattern_16","Pattern_17","Pattern_18","Pattern_19",
                "Pattern_20","Pattern_21")


#-- FIG 4A
m <- multiv_coxph(covariates, no_norm) 
#-- Order factor levels alphabetically
m$pattern <- as.character(rownames(m))
m$pattern <- factor(m$pattern, levels=unique(m$pattern))
m$sig <- m$`Pr(>|z|)` < 0.05
m$se <- m$`se(coef)`


tiff("tcga_legacy_os_by_pattern_scaledpoints.tiff", units="in", width=8, height=3, res=300)
ggplot(m, aes(x = pattern, y = coef)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = pattern, 
                 y = coef, color = sig, size = -log10(`Pr(>|z|)`)), position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  xlab("Pattern") + ylab("Regression Coefficient")
dev.off()

#-- Model cox survival for each cancer type by pattern
#-- Age as covariates with each pattern

cancer_types <- unique(no_norm$cancer)
cancer_multiv_list <- list() 
for (cancer in cancer_types){
  cancer_data <- no_norm[no_norm$cancer == cancer,]
  print(cancer)
  cancer_multiv_list[[cancer]] <- multiv_coxph(covariates, cancer_data) 
}


#-- FIG 4B
#--  pull out specific pattern across each cancer
my_pattern <- lapply(cancer_multiv_list,
       function(x){ 
         pattern_row <- x[7,]
         return(pattern_row)
       })
res <- do.call(rbind.data.frame, my_pattern)
res$sig <- res$`Pr(>|z|)` <= 0.05
res$se <- res$`se(coef)`
res$cancer <- rownames(res)

tiff("tcga_legacy_os_pattern7_by_cancer_scaledpoints.tiff", units="in", width=8, height=3, res=300)
ggplot(res, aes(x = cancer, y = coef)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = cancer, 
                 y = coef, color = sig, size = -log10(`Pr(>|z|)`)), position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  xlab("Cancer") + ylab("Regression Coefficient")
dev.off()




#-- Model cox survival for each cancer type and category by pattern
#-- Age as covariates with each pattern
primary <- no_norm[no_norm$category == "Primary",]
cancer_types <- unique(primary$cancer)

primary_multiv_list <- list() 
for (cancer in cancer_types){
  cancer_data <- primary[primary$cancer == cancer,]
  print(cancer)
  primary_multiv_list[[cancer]] <- multiv_coxph(covariates, cancer_data) 
}

#--  pull out specific pattern across each cancer
my_pattern <- lapply(primary_multiv_list,
                     function(x){ 
                       pattern_row <- x[7,]
                       return(pattern_row)
                     })
res <- do.call(rbind.data.frame, my_pattern)
res$sig <- res$`Pr(>|z|)` <= 0.05
res$se <- res$`se(coef)`
res$cancer <- rownames(res)

tiff("tcga_legacy_os_pattern7_primary_by_cancer_scaledpoints.tiff", units="in", width=8, height=3, res=300)
ggplot(res, aes(x = cancer, y = coef)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = cancer, 
                 y = coef, color = sig, size = -log10(`Pr(>|z|)`)), position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  xlab("Cancer") + ylab("Regression Coefficient")
dev.off()

#-- For melanoma only 

met <- no_norm[no_norm$category == "Metastatic",]
met_skcm <- met[met$cancer == "TCGA-SKCM",]
skcm <- multiv_coxph(covariates, met_skcm)

#-- Order factor levels alphabetically
skcm$pattern <- as.character(rownames(skcm))
skcm$pattern <- factor(skcm$pattern, levels=unique(skcm$pattern))
skcm$sig <- skcm$`Pr(>|z|)` <= 0.05
skcm$se <- skcm$`se(coef)`

tiff("tcga_legacy_os_pattern7_skcm_metastatic_scaledpoints.tiff", units="in", width=8, height=3, res=300)
ggplot(skcm, aes(x = pattern, y = coef)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = pattern, 
                 y = coef, color = sig, size = -log10(`Pr(>|z|)`)), position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  xlab("Pattern") + ylab("Regression Coefficient")
dev.off()

primary <- no_norm[no_norm$category == "Primary",]
primary_skcm <- primary[primary$cancer == "TCGA-SKCM",]
skcm <- multiv_coxph(covariates, primary_skcm)

#-- Order factor levels alphabetically
skcm$pattern <- as.character(rownames(skcm))
skcm$pattern <- factor(skcm$pattern, levels=unique(skcm$pattern))
skcm$sig <- skcm$`Pr(>|z|)` <= 0.05
skcm$se <- skcm$`se(coef)`

tiff("tcga_legacy_os_pattern7_skcm_primary_scaledpoints.tiff", units="in", width=8, height=3, res=300)
ggplot(skcm, aes(x = pattern, y = coef)) +
  scale_color_manual(values = c("black","red")) +
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
  geom_point(aes(x = pattern, 
                 y = coef, color = sig,  size = -log10(`Pr(>|z|)`)), position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=coef-se, ymax=coef+se), width=.2,
                position=position_dodge(0.05)) +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) +
  xlab("Pattern") + ylab("Regression Coefficient")
dev.off()

#-- Plot for hazard ratios
fit <- coxph(Surv(OS.time, OS)~Age+Pattern_7, data = met_skcm)
ggforest(fit, data = met_skcm)



