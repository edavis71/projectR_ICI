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
#library(projectR)
#library(CePa)
#library(org.Hs.eg.db)
#library(biomaRt)
#library(gplots)
library(reshape2)
library(ggplot2)
#library(CoGAPS)
#library(data.table)
library(ComplexHeatmap)
library(viridis)
#library(GEOquery)
library(RColorBrewer)
#library(ROCR)
library(dplyr)
#library(ggalluvial)
library(TCGAbiolinks)
#library(SummarizedExperiment)
library(readxl)
library(car)
library(caret)
library(randomForest)
library(survminer)
library(RTCGA.clinical)
library(survival)
library(MLeval)

# set working directory
setwd('/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI')

# get valid project IDs
project_IDs <- TCGAbiolinks:::getGDCprojects()$project_id
# get "TCGA" named projects
cancer_types <- project_IDs[grep("TCGA", TCGAbiolinks:::getGDCprojects()$project_id)]

###################################################
### read in data from rds files I made (per cancer type if needed)
###################################################
expTCGA <- readRDS("TCGA_RDS/TCGA.legacy.expression.rds")
metaTCGA <- readRDS("TCGA_RDS/TCGA.legacy.meta.rds")

#############################################################################
# Load CoGAPS results object from Gubin et al. 
# and run transfer learning with projectR
#############################################################################
## load CoGAPS result object
# load("myeloid_run7_result3.RData")
# 
# # remove gene duplicates from TCGA
# expTCGA <- expTCGA[match(unique(rownames(expTCGA)), rownames(expTCGA)),]
# # remove gene duplicates from CoGAPs feature loadings
# gapsResult@featureLoadings <- gapsResult@featureLoadings[
#   match(unique(rownames(gapsResult@featureLoadings)),
#         rownames(gapsResult@featureLoadings)),]
# # restrict analysis to the genes that interset between TCGA
# # and gubin CoGAPs run
# # first remove duplicate gene names from both
# symbols <- intersect(rownames(expTCGA),
#               rownames(gapsResult@featureLoadings))
# # subset TCGA and Gubin to intersecting genes
# expTCGA <- expTCGA[rownames(expTCGA) %in% symbols,]
# gapsResult@featureLoadings <- gapsResult@featureLoadings[rownames(gapsResult@featureLoadings) %in% symbols,]
# # run TL
# gaps2TCGA <- projectR(data = expTCGA, 
#                       loadings = gapsResult,
#                       full = TRUE, 
#                       dataNames = symbols)
# ## get projected matrix out of new TL object
# TL.proj <- gaps2TCGA[["projection"]]
# ## set any NAs to zero
# TL.proj[which(is.na(TL.proj))] <- 0
# ## confirm NAs were removed and no other atypical value types (Inf, Nan, etc)
# stopifnot(!is.na(TL.proj))
# stopifnot(!is.infinite(TL.proj))
# stopifnot(!is.nan(TL.proj))

## save TL results
#saveRDS(TL.proj, file = "TL.TCGA.rds")
# read in if starting analysis from here and don't want to rerun steps above
TL.proj <- readRDS("TL.TCGA.rds")

############################################################################
# Process TL output and meta data
#############################################################################
#transpose TL.proj
TL.proj <- as.data.frame(t(TL.proj))

# add relevant meta data to TL.prov df that contains projected patterns
mtch <- match(rownames(TL.proj), metaTCGA$barcode) # make sure dfs in same order
metaTCGA <- metaTCGA[mtch,]
TL.proj$cancertype <- as.factor(metaTCGA$cancertype)

## add survival measures from Liu et al 2018 https://www.sciencedirect.com/science/article/pii/S0092867418302290#app2
survival_data <- read_excel("1-s2.0-S0092867418302290-mmc1.xlsx")

# add age and survival data to TL.proj df
TL.proj$barcode <- sapply(strsplit(rownames(TL.proj), "-", fixed = T),
                          function(x) paste(x[1:3], collapse = "-"))
mtch <- match(TL.proj$barcode, survival_data$bcr_patient_barcode)
# add age at diagnosis
TL.proj$age <- survival_data$age_at_initial_pathologic_diagnosis[mtch]
# add survival measures
TL.proj$OS.time <- survival_data$OS.time[mtch]
TL.proj$DSS.time <- survival_data$DSS.time[mtch]
TL.proj$DFI.time <- survival_data$DFI.time[mtch]
TL.proj$PFI.time <- survival_data$PFI.time[mtch]
# add race and gender
TL.proj$gender <- survival_data$gender[mtch]
TL.proj$race <- survival_data$race[mtch]
# add name
TL.proj$name <- rownames(TL.proj)

############################################################################
# Run Models And Graph Results
#############################################################################
# center and scale the patterns
TL.proj[1:21] <- scale(TL.proj[1:21], center = TRUE, scale = TRUE)
TL.proj$cancertype <- as.factor(TL.proj$cancertype)
# melt the data
# melt dataset
TL.proj.m <- reshape2::melt(TL.proj,
                            id.vars = c("OS.time", "DSS.time", "DFI.time", "PFI.time",
                                        "gender","race", "cancertype", "age", "name"),
                            na.rm = FALSE, value.name = "value")

# run a linear model with y = survival,
# x = cancer types + every CoGAPs pattern (No Age)
lm.all <- lm(TL.proj$OS.time~TL.proj$cancertype+TL.proj$Pattern_1+
                  TL.proj$Pattern_2+TL.proj$Pattern_3+TL.proj$Pattern_4+
                  TL.proj$Pattern_5+TL.proj$Pattern_6+TL.proj$Pattern_7+
                  TL.proj$Pattern_8+TL.proj$Pattern_9+TL.proj$Pattern_10+
                  TL.proj$Pattern_11+TL.proj$Pattern_12+TL.proj$Pattern_13+
                  TL.proj$Pattern_14+TL.proj$Pattern_15+TL.proj$Pattern_16+
                  TL.proj$Pattern_17+TL.proj$Pattern_18+TL.proj$Pattern_19+
                  TL.proj$Pattern_20+TL.proj$Pattern_21)

# make a plot of normalized coefficients with CI
# first get coefficnents and CIs
lm.all.coef <- summary(lm.all)$coefficients
xvals <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),1]
CIlow <- xvals - (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
CIhigh <- xvals + (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
pval <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),4]
yvals <- gsub("TL.proj$", "", names(xvals), fixed = T)
dat.df <- data.frame(yvals = factor(yvals, levels = yvals),
                     xvals = xvals,
                     CIlow = CIlow,
                     CIhigh = CIhigh,
                     pval = pval)
# make plot using ggplot
plotcols <- rep("black", nrow(dat.df))
plotcols[dat.df$pval <= 0.05] <- "red"
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = plotcols) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
            #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Standardized Coefficient")
  #annotate(geom = "text", y =1.1, x = log10(1.5), 
           #label = "Model p < 0.001\nPseudo R^2 = 0.10", size = 3.5, hjust = 0) + 
  #ggtitle("Feeding method and risk of obesity in cats")

pdf("ICI_projectR.TCGA.lm.all.no_age.pdf")
p
dev.off()

## linear model WITH age
lm.all <- lm(TL.proj$OS.time~TL.proj$cancertype+TL.proj$Pattern_1+
               TL.proj$Pattern_2+TL.proj$Pattern_3+TL.proj$Pattern_4+
               TL.proj$Pattern_5+TL.proj$Pattern_6+TL.proj$Pattern_7+
               TL.proj$Pattern_8+TL.proj$Pattern_9+TL.proj$Pattern_10+
               TL.proj$Pattern_11+TL.proj$Pattern_12+TL.proj$Pattern_13+
               TL.proj$Pattern_14+TL.proj$Pattern_15+TL.proj$Pattern_16+
               TL.proj$Pattern_17+TL.proj$Pattern_18+TL.proj$Pattern_19+
               TL.proj$Pattern_20+TL.proj$Pattern_21+TL.proj$age)

# make a plot of normalized coefficients with CI
# first get coefficnents and CIs
lm.all.coef <- summary(lm.all)$coefficients
xvals <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),1]
CIlow <- xvals - (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
CIhigh <- xvals + (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
pval <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),4]
yvals <- gsub("TL.proj$", "", names(xvals), fixed = T)
dat.df <- data.frame(yvals = factor(yvals, levels = yvals),
                     xvals = xvals,
                     CIlow = CIlow,
                     CIhigh = CIhigh,
                     pval = pval)
# make plot using ggplot
plotcols <- rep("black", nrow(dat.df))
plotcols[dat.df$pval <= 0.05] <- "red"
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = plotcols) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Standardized Coefficient")

pdf("ICI_projectR.TCGA.lm.all.age.pdf")
p
dev.off()

# association between patterns and age
lm.all <- lm(TL.proj$age~TL.proj$cancertype+TL.proj$Pattern_1+
               TL.proj$Pattern_2+TL.proj$Pattern_3+TL.proj$Pattern_4+
               TL.proj$Pattern_5+TL.proj$Pattern_6+TL.proj$Pattern_7+
               TL.proj$Pattern_8+TL.proj$Pattern_9+TL.proj$Pattern_10+
               TL.proj$Pattern_11+TL.proj$Pattern_12+TL.proj$Pattern_13+
               TL.proj$Pattern_14+TL.proj$Pattern_15+TL.proj$Pattern_16+
               TL.proj$Pattern_17+TL.proj$Pattern_18+TL.proj$Pattern_19+
               TL.proj$Pattern_20+TL.proj$Pattern_21)

# make a plot of normalized coefficients with CI
# first get coefficnents and CIs
lm.all.coef <- summary(lm.all)$coefficients
xvals <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),1]
CIlow <- xvals - (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
CIhigh <- xvals + (1.96*lm.all.coef[grep("Pattern", rownames(lm.all.coef)),2])
pval <- lm.all.coef[grep("Pattern", rownames(lm.all.coef)),4]
yvals <- gsub("TL.proj$", "", names(xvals), fixed = T)
dat.df <- data.frame(yvals = factor(yvals, levels = yvals),
                     xvals = xvals,
                     CIlow = CIlow,
                     CIhigh = CIhigh,
                     pval = pval)
# make plot using ggplot
plotcols <- rep("black", nrow(dat.df))
plotcols[dat.df$pval <= 0.05] <- "red"
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = plotcols) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Standardized Coefficient")

pdf("ICI_projectR.TCGA.lm.age_vs_patterns.pdf")
p
dev.off()

# run linear models per cancer type
ct.names <- c()
ct.out <- c()
for (ct in levels(TL.proj$cancertype)){
  TL.proj.ct <- subset(TL.proj, cancertype == ct)
  lm.ct <- lm(TL.proj.ct$OS.time~TL.proj.ct$Pattern_1+
                TL.proj.ct$Pattern_2+TL.proj.ct$Pattern_3+
                TL.proj.ct$Pattern_4+TL.proj.ct$Pattern_5+
                TL.proj.ct$Pattern_6+TL.proj.ct$Pattern_7+
                TL.proj.ct$Pattern_8+TL.proj.ct$Pattern_9+
                TL.proj.ct$Pattern_10+TL.proj.ct$Pattern_11+
                TL.proj.ct$Pattern_12+TL.proj.ct$Pattern_13+
                TL.proj.ct$Pattern_14+TL.proj.ct$Pattern_15+
                TL.proj.ct$Pattern_16+TL.proj.ct$Pattern_17+
                TL.proj.ct$Pattern_18+TL.proj.ct$Pattern_19+
                TL.proj.ct$Pattern_20+TL.proj.ct$Pattern_21+
                TL.proj.ct$age)
  #print_output <- paste0(ct,": ",rownames(summary(lm.ct)$coefficients)[8]," ",paste0(summary(lm.ct)$coefficients[8,],collapse = " "))
  ct.names <- append(ct.names, ct)
  ct.out <- rbind(ct.out, summary(lm.ct)$coefficients[8,])
  #print(print_output)
}
# convert output to df, define CI, and plot
ct.df <- as.data.frame(ct.out)
ct.df$Names <- ct.names
ct.df$CIlow <- ct.df$Estimate - (1.96*ct.df$`Std. Error`)
ct.df$CIhigh <- ct.df$Estimate + (1.96*ct.df$`Std. Error`)
colnames(ct.df)[4] <- "pval"
# plot results
plotcols <- rep("black", nrow(ct.df))
plotcols[ct.df$pval <= 0.05] <- "red"
p <- ggplot(ct.df, aes(x = Estimate, y = Names, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = plotcols) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Standardized Coefficient")

# save at size height = 14, width = 7
pdf("ICI_projectR.TCGA.P7_all.pdf", height = 14, width = 7)
p
dev.off()

# explore patterns within SKCM in great depth
TL.proj.SKCM <- subset(TL.proj, cancertype == "TCGA-SKCM")
lm.SKCM <- lm(TL.proj.SKCM$OS.time~TL.proj.SKCM$Pattern_1+
                    TL.proj.SKCM$Pattern_2+TL.proj.SKCM$Pattern_3+TL.proj.SKCM$Pattern_4+
                    TL.proj.SKCM$Pattern_5+TL.proj.SKCM$Pattern_6+TL.proj.SKCM$Pattern_7+
                    TL.proj.SKCM$Pattern_8+TL.proj.SKCM$Pattern_9+TL.proj.SKCM$Pattern_10+
                    TL.proj.SKCM$Pattern_11+TL.proj.SKCM$Pattern_12+TL.proj.SKCM$Pattern_13+
                    TL.proj.SKCM$Pattern_14+TL.proj.SKCM$Pattern_15+TL.proj.SKCM$Pattern_16+
                    TL.proj.SKCM$Pattern_17+TL.proj.SKCM$Pattern_18+TL.proj.SKCM$Pattern_19+
                    TL.proj.SKCM$Pattern_20+TL.proj.SKCM$Pattern_21+TL.proj.SKCM$age)

# plot coefficients from SKCM model
# first get coefficnents and CIs
lm.SKCM.coef <- summary(lm.SKCM)$coefficients
xvals <- lm.SKCM.coef[grep("Pattern", rownames(lm.SKCM.coef)),1]
CIlow <- xvals - (1.96*lm.SKCM.coef[grep("Pattern", rownames(lm.SKCM.coef)),2])
CIhigh <- xvals + (1.96*lm.SKCM.coef[grep("Pattern", rownames(lm.SKCM.coef)),2])
pval <- lm.SKCM.coef[grep("Pattern", rownames(lm.SKCM.coef)),4]
yvals <- gsub("TL.proj.SKCM$", "", names(xvals), fixed = T)
yvals_order <- yvals
dat.df <- data.frame(yvals = factor(yvals, levels = yvals_order),
                     xvals = xvals,
                     CIlow = CIlow,
                     CIhigh = CIhigh,
                     pval = pval)

# make plot using ggplot
plotcols <- rep("black", nrow(dat.df))
plotcols[dat.df$pval <= 0.05] <- "red"
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = plotcols) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Standardized Coefficient")

# save at size height = 15, width = 5
pdf("ICI_projectR.TCGA.lm.SKCM.h15.w5.pdf", height = 15, width = 5)
p
dev.off()
# save at size height = 15, width = 10
pdf("ICI_projectR.TCGA.lm.SKCM.h15.w10.pdf", height = 15, width = 10)
p
dev.off()
# save at size height = 15, width = 8
pdf("ICI_projectR.TCGA.lm.SKCM.h15.w8.pdf", height = 15, width = 8)
p
dev.off()
# save at size height = 12, width = 8
pdf("ICI_projectR.TCGA.lm.SKCM.h12.w8.pdf", height = 12, width = 8)
p
dev.off()
# save at size height = 15, width = 6
pdf("ICI_projectR.TCGA.lm.SKCM.h15.w6.pdf", height = 15, width = 6)
p
dev.off()
# save at size height = 12, width = 6
pdf("ICI_projectR.TCGA.lm.SKCM.h12.w6.pdf", height = 12, width = 6)
p
dev.off()

# Correlation betwen B7 and patterns in all cancers and SKCM specifically
expTCGA.t <- t(expTCGA)
expTCGA.t.B7 <- expTCGA.t[,which(colnames(expTCGA.t) == "CD80" | colnames(expTCGA.t) == "CD86")]
mtch <- match(rownames(TL.proj), rownames(expTCGA.t.B7))
TL.proj.B7 <- cbind(TL.proj, expTCGA.t.B7[mtch,])

# all TCGA samples and B7.1
patterns.B7.1.all <- apply(TL.proj.B7[,c(1:21)], 2,
              function(a){
                res <- cor.test(a, TL.proj.B7[,32])
                return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.1.df <- as.data.frame(t(patterns.B7.1.all))
colnames(patterns.B7.1.df) <- c("Estimate","CIlow",
                                 "CIhigh", "pval")
patterns.B7.1.df$pattern <- rownames(patterns.B7.1.df)
patterns.B7.1.df$plotcol <- "black"
patterns.B7.1.df$plotcol[patterns.B7.1.df$pval <= 0.05] <- "red"
# write data to outfile
write.csv(patterns.B7.1.df, file = "patterns.B7.1.csv",
          quote = F, row.names = F)
# plot
# make plot using ggplot
p <- ggplot(patterns.B7.1.df, aes(x = Estimate, y = pattern,
                                  size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "grey") +
  geom_point(color = patterns.B7.1.df$plotcol) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Pearson Correlation Coefficient")

pdf("patterns.B7.1.Pearson.pdf", height = 12, width = 6)
p
dev.off()

# all TCGA samples and B7.2

patterns.B7.2.all <- apply(TL.proj.B7[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7[,33])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.2.df <- as.data.frame(t(patterns.B7.2.all))
colnames(patterns.B7.2.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
# if pval is zero, set to min for graphing
patterns.B7.2.df$pval[patterns.B7.2.df$pval == 0] <- min(patterns.B7.2.df$pval[patterns.B7.2.df$pval != 0])

patterns.B7.2.df$pattern <- rownames(patterns.B7.2.df)
patterns.B7.2.df$plotcol <- "black"
patterns.B7.2.df$plotcol[patterns.B7.2.df$pval <= 0.05] <- "red"
#View(patterns.B7.2.df)
# write data to outfile
write.csv(patterns.B7.2.df, file = "patterns.B7.2.csv",
          quote = F, row.names = F)
# plot
# make plot using ggplot
p <- ggplot(patterns.B7.2.df, aes(x = Estimate, y = pattern,
                                  size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "grey") +
  geom_point(color = patterns.B7.2.df$plotcol) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Pearson Correlation Coefficient")

pdf("patterns.B7.2.Pearson.pdf", height = 12, width = 6)
p
dev.off()

# now repeat but just with SKCM
TL.proj.B7.SKCM <- TL.proj.B7[TL.proj.B7$cancertype == "TCGA-SKCM",]

# all TCGA samples and B7.1
patterns.B7.1.SKCM <- apply(TL.proj.B7.SKCM[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7.SKCM[,32])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.1.SKCM.df <- as.data.frame(t(patterns.B7.1.SKCM))
colnames(patterns.B7.1.SKCM.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
patterns.B7.1.SKCM.df$pattern <- rownames(patterns.B7.1.SKCM.df)
patterns.B7.1.SKCM.df$plotcol <- "black"
patterns.B7.1.SKCM.df$plotcol[patterns.B7.1.SKCM.df$pval <= 0.05] <- "red"
# write data to outfile
write.csv(patterns.B7.1.SKCM.df, file = "patterns.B7.1.SKCM.csv",
          quote = F, row.names = F)
# plot
# make plot using ggplot
p <- ggplot(patterns.B7.1.SKCM.df, aes(x = Estimate, y = pattern,
                                  size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "grey") +
  geom_point(color = patterns.B7.1.SKCM.df$plotcol) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Pearson Correlation Coefficient")

pdf("patterns.B7.1.SKCM.Pearson.pdf", height = 12, width = 6)
p
dev.off()

# all TCGA samples and B7.2

patterns.B7.2.SKCM <- apply(TL.proj.B7.SKCM[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7.SKCM[,33])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.2.SKCM.df <- as.data.frame(t(patterns.B7.2.SKCM))
colnames(patterns.B7.2.SKCM.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
# if pval is zero, set to min for graphing
patterns.B7.2.SKCM.df$pval[patterns.B7.2.SKCM.df$pval == 0] <- min(patterns.B7.2.SKCM.df$pval[patterns.B7.2.SKCM.df$pval != 0])

patterns.B7.2.SKCM.df$pattern <- rownames(patterns.B7.2.SKCM.df)
patterns.B7.2.SKCM.df$plotcol <- "black"
patterns.B7.2.SKCM.df$plotcol[patterns.B7.2.SKCM.df$pval <= 0.05] <- "red"
#View(patterns.B7.2.df)
# write data to outfile
write.csv(patterns.B7.2.SKCM.df, file = "patterns.B7.2.SKCM.csv",
          quote = F, row.names = F)
# plot
# make plot using ggplot
p <- ggplot(patterns.B7.2.SKCM.df, aes(x = Estimate, y = pattern,
                                  size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "grey") +
  geom_point(color = patterns.B7.2.SKCM.df$plotcol) +
  #coord_trans(x = scales:::exp_trans(10)) +
  #scale_x_continuous(breaks = log10(seq(0.1, 2.5, 0.1)), labels = seq(0.1, 2.5, 0.1),
  #         limits = log10(c(0.09,2.5))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  ylab("") +
  xlab("Pearson Correlation Coefficient")

pdf("patterns.B7.2.SKCM.Pearson.pdf", height = 12, width = 6)
p
dev.off()

# comment out barplots
# pdf("patterns.B7.1.all.pdf")
# barplot(patterns.B7.1.all[,1], las = 2,
#         main = "patterns.B7.1.all", ylab = "Pearson Correlation")
# dev.off()
# 
# pdf("patterns.B7.2.all.pdf")
# barplot(patterns.B7.2.all[,1], las = 2,
#         main = "patterns.B7.2.all", ylab = "Pearson Correlation")
# dev.off()
# 
# pdf("patterns.B7.1.SKCM.pdf")
# barplot(patterns.B7.1.SKCM[,1], las = 2,
#         main = "patterns.B7.1.SKCM", ylab = "Pearson Correlation")
# dev.off()
# 
# pdf("patterns.B7.2.SKCM.pdf")
# barplot(patterns.B7.2.SKCM[,1], las = 2,
#         main = "patterns.B7.2.SKCM", ylab = "Pearson Correlation")
# dev.off()

# Note: The analysis below (using estimates from xCell and CIBERSORTx)
# did not identify any interesting relationships
# One potential reason for this is that the estimates from CIBERSORTx
# are too noisy, especially when using the standard 22 immune cell
# reference set used in this CIBERSORTx analysis

# Compare activated and resting NK cell proportions, as estimated by CIBERSORTx,
# with TL pattern weights
# TCGA.cbr <- read.table("CIBERSORTx_SKCM_Adjusted.txt", sep = "\t", header = T)
# #TCGA.cbr <- subset(TCGA.cbr, Correlation >= 0.7)
# # save rownames
# TCGA.cbr.rownames <- TCGA.cbr[,1]
# # select cols with cell proportion estimates
# TCGA.cbr <- TCGA.cbr[,2:23]
# # assign rownames
# rownames(TCGA.cbr) <- TCGA.cbr.rownames
# # # match tcga and cibersort rows
# mtch <- match(rownames(TL.proj), rownames(TCGA.cbr))
# TL.proj.cbr <- TL.proj[which(!is.na(mtch)),]
# TL.proj.cbr$NKa <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$NK.cells.activated)
# TL.proj.cbr$NKr <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$NK.cells.resting)
# TL.proj.cbr$Tregs <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$T.cells.regulatory..Tregs.)
# TL.proj.cbr$T8 <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$T.cells.CD8)
# TL.proj.cbr$MacM1 <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$Macrophages.M1)
# TL.proj.cbr$MacM2 <- as.numeric(TCGA.cbr[mtch[which(!is.na(mtch))],]$Macrophages.M2)
# 
# View(TL.proj.cbr)
# View(head(expTCGA))
# CTLA4.exp <- expTCGA[which(rownames(expTCGA) == "CTLA4"),]
# mtch <- match(rownames(TL.proj.cbr), names(CTLA4.exp))
# TL.proj.cbr$CTLA4 <- CTLA4.exp[mtch]
# TL.proj.cbr.sub <- subset(TL.proj.cbr, NKa >= 0.01)
# # 
# cor(TL.proj.cbr$NKa,TL.proj.cbr$CTLA4, method = "spearman")
# cor(TL.proj.cbr$Pattern_7,TL.proj.cbr$CTLA4)
# cor(TL.proj.cbr.sub$NKr,TL.proj.cbr.sub$CTLA4)

# 
# ggplot(TL.proj.cbr.sub,
#        aes(x = Pattern_7, y = NKr)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "red") +
#   theme_classic()
  #facet_wrap(~cancertype, scales = "free")

# compare xCell proportions and patterns using published xCell estimates for TCGA
# read TCGA data
# TCGA.xCell <- read.table("xCell_TCGA_RSEM.txt", sep = "\t", header = T)
# # save sample names
# TCGA.xCell.colnames <- TCGA.xCell$X
# # remove col X before normalization (contains sample names)
# TCGA.xCell <- TCGA.xCell[,2:ncol(TCGA.xCell)]
# # normalize values per sample
# TCGA.xCell <- apply(TCGA.xCell, 2, function(x) as.numeric(as.character(x))/sum(as.numeric(as.character(x)),na.rm = T))
# # traspose xCell dfs
# TCGA.xCell <- t(TCGA.xCell)
# TCGA.xCell <- as.data.frame(TCGA.xCell)
# View(TCGA.xCell)
# # only keep TCGA samples
# TCGA.xCell <- TCGA.xCell[grepl("TCGA", rownames(TCGA.xCell)),]
# # reassign sample names as rownames
# rownames(TCGA.xCell) <- gsub(".", "-",
#                              substring(rownames(TCGA.xCell), 1, 12),
#                              fixed = T)
# # set col names
# colnames(TCGA.xCell) <- TCGA.xCell.colnames
# # match tcga and xcell rows
# mtch <- match(TL.proj$barcode, rownames(TCGA.xCell))
# TL.proj$NK <- TCGA.xCell[mtch,]$`NK cells`
# TL.proj$NKT <- TCGA.xCell[mtch,]$NKT
# TL.proj$Tregs <- TCGA.xCell[mtch,]$Tregs
# TL.proj$Th1 <- TCGA.xCell[mtch,]$`Th1 cells`
# TL.proj.dis <- subset(TL.proj, cancertype == "TCGA-KIRC")
# ggplot(TL.proj,
#         aes(x = Pattern_7, y = NKT)) +
#         geom_point() +
#         geom_smooth(method = "lm", col = "red") +
#         theme_classic() +
#         facet_wrap(~cancertype, scales = "free")
# 
# # compare CIBERSORTx proportions and patterns
# # get cibersort results for TCGA from supp file from Thorsson et al. 2018
# TCGAthorsson <- as.data.frame(read_excel("1-s2.0-S1074761318301213-mmc2.xlsx"))
# # subset and transpose
# TCGAcbr <- TCGAthorsson[,37:58]
# # assign TCGA ID as col name
# rownames(TCGAcbr) <- TCGAthorsson$`TCGA Participant Barcode`
# # recode colnames
# colnames(TCGAcbr) <- sapply(sapply(colnames(TCGAcbr), function(x) strsplit(x, "...", fixed = T)), "[", 1)
# View(rownames(TCGAcbr))
# # match tcga and cibersort rows
# mtch <- match(TL.proj$barcode, rownames(TCGAcbr))
# TL.proj$NKa <- as.numeric(TCGAcbr[mtch,]$`NK Cells Activated`)
# TL.proj$NKr <- as.numeric(TCGAcbr[mtch,]$`NK Cells Resting`)
# TL.proj$Tregs <- as.numeric(TCGAcbr[mtch,]$`T Cells Regulatory Tregs`)
# TL.proj$T8 <- as.numeric(TCGAcbr[mtch,]$`T Cells CD8`)
# TL.proj$MacM2 <- as.numeric(TCGAcbr[mtch,]$`Macrophages M2`)
# TL.proj.dis <- subset(TL.proj, cancertype == "TCGA-KIRC")
# ggplot(TL.proj,
#        aes(x = Pattern_7, y = NKr)) +
#   geom_point() +
#   geom_smooth(method = "lm", col = "red") +
#   theme_classic() +
#   facet_wrap(~cancertype, scales = "free")

  
# explore kaplan meyer survival

# SKCM and pattern 7
# Pat <- subset(TL.proj, cancertype == "TCGA-SKCM")$Pattern_7
# OS <- subset(TL.proj, cancertype == "TCGA-SKCM")$OS.time
# group <- rep(NA, length(OS))
# event <- rep(1, length(OS))
# TL.proj.km <- data.frame(OS = OS, event = event, group = group)
# cutoffs <- quantile(Pat, probs = c(0.35, 0.65))
# # stratify based on pattern 7
# TL.proj.km$group[Pat <= cutoffs[1]] <- "low"
# TL.proj.km$group[Pat >= cutoffs[2]] <- "high"
# #TL.proj.km <- TL.proj.km[!is.na(TL.proj.km$group),]
# #TL.proj.km$group <- as.factor(TL.proj.km$group)
# # calc survival model fit
# KMfit <- survfit(Surv(OS, event) ~ group,
#                data = TL.proj.km)
# # plot with survminer
# ggsurvplot(KMfit, data = TL.proj.km, risk.table = TRUE, pval = T)

# Explore prediction of overall survival using non-linear predictors (random
# forest, MARS, etc)

# Train random forest models for the prediction of overall survival from all cancer types, patterns, and age
# y <- TL.proj$OS.time
# x <- TL.proj[,c(1:22,24)]
# rf.in <- cbind(y,x)
# # scale age - you already scaled the patterns
# rf.in$age <- scale(rf.in$age, center = T, scale = T)
# # remove NAs
# rf.in <- rf.in[-which(apply(rf.in, 1, function(a) sum(is.na(a)) > 0 )),]
# 
# model_rf = train(y ~ ., data=rf.in, method='rf')
# model_rf
# plot(model_rf)
# varimp_mars <- varImp(model_rf)
# plot(varimp_mars, main="Variable Importance with RF")
# 
# testData2 <- predict(preProcess_range_model, testData)
# predicted <- predict(model_rf, testData2)
# 
# plot(testData$Sensitivity, predicted)
# TL.proj.SKCM <- subset(TL.proj, cancertype == "TCGA-SKCM")
# # Train random forest model to predict survival using patterns + age in melanoma (SKCM) tumors
# y <- TL.proj.SKCM$OS.time
# x <- TL.proj.SKCM[,c(1:21,24)]
# rf.in <- cbind(y,x)
# # scale age - you already scaled the patterns
# rf.in$age <- scale(rf.in$age, center = T, scale = T)
# # remove NAs
# rf.in <- rf.in[-which(apply(rf.in, 1, function(a) sum(is.na(a)) > 0 )),]
# model_rf = train(y ~ ., data=rf.in, method='lm')
# model_rf$results
# plot(model_rf)
# varimp_mars <- varImp(model_rf)
# plot(varimp_mars, main="Variable Importance with RF")
# 
# model_rf$finalModel
# 
# # testData2 <- predict(preProcess_range_model, testData)
# predicted <- predict(model_rf, rf.in)
# cor(rf.in$y, predicted)
# 
# myTrainingControl <- trainControl(method = "repeatedcv", 
#                                   number = 10, 
#                                   savePredictions = TRUE,
#                                   repeats = 5,
#                                   #classProbs = TRUE,
#                                   verboseIter = TRUE)
# 
# randomForestFit = train(x = rf.in[,2:23], 
#                         y = rf.in$y, 
#                         method = "rf", 
#                         trControl = myTrainingControl) 
#                         #preProcess = c("center","scale"))
# 
# pdf("ICI_projectR.TCGA.rf.SKCM.varImp.pdf")
# plot(varImp(randomForestFit))
# dev.off()
# 
# randomForestFit$finalModel
# 
# eval.rf <- evalm(randomForestFit)
# 
# ## get roc curve plotted in ggplot2
# 
# x$roc

## get AUC and other metrics

#x$stdres


# # explore clustering and heatmap visualizations - this was unremarkable
# 
# # see if cancers cluster by patterns
# col_vector = color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
# anno_cols <- sample(col_vector, length(levels(TL.proj$cancertype)))
# pdf("col_test.pdf")
# plot(1:length(anno_cols), rep(10,length(anno_cols)), pch = 19, cex = 2,
#      col = anno_cols)
# dev.off()
# anno_list <- list()
# for (j in 1:length(levels(TL.proj$cancertype))){
#   anno_list[[levels(TL.proj$cancertype)[j]]] <- anno_cols[j]
# }
# row_ha = rowAnnotation("Cancer Types" = TL.proj$cancertype,
#                        col = list("Cancer Types" = unlist(anno_list)))
#                          
#                        #   anno_block(gp = gpar(fill = rainbow(length(levels(TL.proj$cancertype)))),
#                        #              labels = levels(TL.proj$cancertype), 
#                        #              labels_gp = gpar(col = "white",
#                        #              fontsize = 10)),
#                        # annotation_name_rot = 90)
# # scale data
# TL.proj.rownames <- rownames(TL.proj)
# TL.proj.scl <- apply(TL.proj[,1:21], 2, scale)
# rownames(TL.proj.scl) <- TL.proj.rownames # reassign rownames
# # open connection to plotting device
# tiff(file = paste("TL.all_CoGAPs_patterns.ComplexHeatmap.TCGA.tiff"),
#      width = 6000, height = 5000, units = "px", res = 800)
# Heatmap(TL.proj.scl, col=inferno(100), name = "mat",
#         clustering_distance_rows = "pearson",
#         column_names_gp = gpar(fontsize = 5),
#         #top_annotation=ha_column, 
#         left_annotation=row_ha,
#         border = TRUE,
#         row_km = 1,
#         show_row_names = FALSE
# )
# dev.off()
# 
# # now see if patterns correlate with time to death
# TL.proj.dtd <- subset(TL.proj, !is.na(paper_Survival))
# # set cols
# anno_list <- list()
# TL.proj.dtd$cancertype <- TL.proj.dtd$cancertype
# for (j in 1:length(levels(TL.proj.dtd$cancertype))){
#   anno_list[[levels(TL.proj.dtd$cancertype)[j]]] <- anno_cols[j]
# }
# 
# row_ha = rowAnnotation("Days to Death" =
#                          anno_lines(TL.proj.dtd$days_to_death,
#                                     axis_param= list(direction = "reverse",
#                                                      labels_rot = 90),
#                                     add_points = F),
#                        "Cancer Types" = TL.proj.dtd$cancertype,
#                        col = list("Cancer Types" = unlist(anno_list)))
# 
#                        # "Immune Score" =
#                        #   anno_lines(TL.proj$ESTIMATE.immune.score,
#                        #              axis_param= list(direction = "reverse",
#                        #                               labels_rot = 90),
#                        #              add_points = F),
# 
# # scale data
# TL.proj.dtd.rownames <- rownames(TL.proj.dtd)
# TL.proj.dtd.scl <- apply(TL.proj.dtd[,1:21], 2, scale)
# rownames(TL.proj.dtd.scl) <- TL.proj.dtd.rownames # reassign rownames
# # open connection to plotting device
# tiff(file = paste("TL.all_CoGAPs_patterns.ComplexHeatmap.TCGA.dtd.tiff"),
#      width = 6000, height = 5000, units = "px", res = 800)
# Heatmap(TL.proj.dtd.scl, col=inferno(100), name = "mat",
#         clustering_distance_rows = "pearson",
#         column_names_gp = gpar(fontsize = 5),
#         #top_annotation=ha_column, 
#         left_annotation=row_ha,
#         border = TRUE,
#         row_km = 1,
#         show_row_names = FALSE
# )
# dev.off()
# 
# # correlation between days to death (dtd) and patterns
# View(TL.proj.scl)
# TL.proj.scl <- as.data.frame(TL.proj.scl)
# TL.proj.scl$immune <- TL.proj$ESTIMATE.immune.score
# TL.proj.scl$days_to_death <- TL.proj$days_to_death
# TL.proj.scl$age <- TL.proj$age_at_index
# TL.proj.scl$name <- rownames(TL.proj.scl)
# TL.proj.scl$cancertype <- TL.proj$cancertype
# TL.proj.scl$cols <- anno_cols[TL.proj$cancertype]
# TL.proj.scl.m <- reshape2::melt(TL.proj.scl,
#       id.vars = c("immune", "days_to_death","age", "cancertype", "cols", "name"),
#       na.rm = FALSE, value.name = "value")
# 
# # explore each pattern
# TL.proj.scl.m$value <- as.numeric(TL.proj.scl.m$value)
# TL.proj.scl.m.pat <- subset(TL.proj.scl.m, variable == "Pattern_4")
# pdf("TL.TCGA.pattern_vs_outcomes.pdf")
# ggplot(TL.proj.scl.m.pat,
#        aes(x = value, y = age, col = "black")) +
#        geom_point() +
#        geom_smooth(method = "lm", se = FALSE) +
#        scale_colour_manual(values = anno_cols) +
#        theme_classic() +
#        facet_wrap(~cancertype, scales = "free")
# dev.off()
# 
