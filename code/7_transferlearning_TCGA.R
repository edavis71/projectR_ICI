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
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(TCGAbiolinks)
library(readxl)
library(car)
library(caret)
library(randomForest)
library(survminer)
library(RTCGA.clinical)
library(survival)
library(MLeval)

# source functions
source("lmPlot.R")

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
load("myeloid_run7_result3.RData")

# remove gene duplicates from TCGA
 expTCGA <- expTCGA[match(unique(rownames(expTCGA)), rownames(expTCGA)),]
# # remove gene duplicates from CoGAPs feature loadings
 gapsResult@featureLoadings <- gapsResult@featureLoadings[
   match(unique(rownames(gapsResult@featureLoadings)),
         rownames(gapsResult@featureLoadings)),]
# restrict analysis to the genes that interset between TCGA
# and gubin CoGAPs run
# first remove duplicate gene names from both
 symbols <- intersect(rownames(expTCGA),
               rownames(gapsResult@featureLoadings))
# subset TCGA and Gubin to intersecting genes
expTCGA <- expTCGA[rownames(expTCGA) %in% symbols,]
gapsResult@featureLoadings <- gapsResult@featureLoadings[rownames(gapsResult@featureLoadings) %in% symbols,]
# # run TL
gaps2TCGA <- projectR(data = expTCGA,
                       loadings = gapsResult,
                       full = TRUE,
                       dataNames = symbols)
# ## get projected matrix out of new TL object
TL.proj <- gaps2TCGA[["projection"]]
# ## set any NAs to zero
TL.proj[which(is.na(TL.proj))] <- 0
# ## confirm NAs were removed and no other atypical value types (Inf, Nan, etc)
stopifnot(!is.na(TL.proj))
stopifnot(!is.infinite(TL.proj))
stopifnot(!is.nan(TL.proj))

## save TL results
saveRDS(TL.proj, file = "TL.TCGA.rds")
# read in if starting analysis from here and don't want to rerun steps above
TL.proj <- readRDS("TL.TCGA.rds")
############################################################################
# Process TL output and meta data
#############################################################################
#transpose TL.proj
TL.proj <- as.data.frame(t(TL.proj))
# filter out adjacent normal samples

# add relevant meta data to TL.prov df that contains projected patterns
mtch <- match(rownames(TL.proj), metaTCGA$barcode) # make sure dfs in same order
metaTCGA <- metaTCGA[mtch,]
# remove NAs if they exits
if (sum(is.na(mtch)) > 0){
  TL.proj <- TL.proj[-which(is.na(mtch)),]
  metaTCGA <- metaTCGA[-which(is.na(mtch)),]
}
TL.proj$cancertype <- as.factor(metaTCGA$cancertype)
TL.proj$definition <- as.factor(metaTCGA$definition)

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
# center and scale the patterns
TL.proj[1:21] <- scale(TL.proj[1:21], center = TRUE, scale = TRUE)

# melt data
# TL.proj.m <- reshape2::melt(TL.proj,
#                             id.vars = c("OS.time", "DSS.time", "DFI.time", "PFI.time",
#                                         "gender","race", "cancertype", "age", "name"),
#                             na.rm = FALSE, value.name = "value")

############################################################################
# Run Models And Graph Results
#############################################################################

# subset samples to run models with only TCGA tumors called:
# Primary solid Tumor, Recurrent Solid Tumor, Metastatic
TL.proj.sub <- subset(TL.proj, definition == "Primary solid Tumor" |
                        definition == "Recurrent Solid Tumor" |
                        definition == "Metastatic")

# run a linear model with y = survival,
# x = cancer types + every CoGAPs pattern (No Age)

mod_cols <- c("cancertype","Pattern_1","Pattern_2","Pattern_3",
              "Pattern_4","Pattern_5","Pattern_6","Pattern_7",
              "Pattern_8","Pattern_9","Pattern_10","Pattern_11",
              "Pattern_12","Pattern_13","Pattern_14","Pattern_15",
              "Pattern_16","Pattern_17","Pattern_18","Pattern_19",
              "Pattern_20","Pattern_21")
                
p <- lmPlot(TL.proj.sub,
            o = "OS.time",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")

# save plot  
pdf("ICI_projectR.TCGA.lm.all.no_age.pdf", height = 12, width = 6)
p
dev.off()

## now run a linear model WITH age as a covariate
mod_cols <- c(mod_cols, "age")
p <- lmPlot(TL.proj.sub,
            o = "OS.time",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")

# save plot
pdf("ICI_projectR.TCGA.lm.all.age.pdf", height = 12, width = 6)
p
dev.off()

# test association between patterns and age
mod_cols <- mod_cols[-length(mod_cols)]
p <- lmPlot(TL.proj.sub,
            o = "age",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")

pdf("ICI_projectR.TCGA.lm.age_vs_patterns.pdf", height = 12, width = 6)
p
dev.off()

# run linear models per cancer type and plot pattern 7 across cancer types
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

# plot pattern 7 across cancer definition types
TL.proj.sub$definition <- factor(TL.proj.sub$definition, levels = unique(as.character(TL.proj.sub$definition)))
# Compare pattern 7 acrosscancer definition types
pdf("ICI_projectR.TCGA.P7.definition.boxplots.pdf")
ggplot(TL.proj.sub,
       aes(x = definition, y = Pattern_7))+
  geom_boxplot() + 
  theme_classic()
dev.off()

# explore patterns within SKCM in great depth
TL.proj.SKCM <- subset(TL.proj.sub, cancertype == "TCGA-SKCM")

# first, plot pattern 7 by definition for only KSCM
pdf("ICI_projectR.TCGA.SKCM.P7.definition.boxplots.pdf")
ggplot(TL.proj.SKCM,
       aes(x = definition, y = Pattern_7))+
  geom_boxplot() + 
  theme_classic()
dev.off()

# now lm model for only SKCM
mod_cols <- c("Pattern_1","Pattern_2","Pattern_3",
              "Pattern_4","Pattern_5","Pattern_6","Pattern_7",
              "Pattern_8","Pattern_9","Pattern_10","Pattern_11",
              "Pattern_12","Pattern_13","Pattern_14","Pattern_15",
              "Pattern_16","Pattern_17","Pattern_18","Pattern_19",
              "Pattern_20","Pattern_21", "age")

p <- lmPlot(TL.proj.SKCM,
            o = "OS.time",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")

# save at size height = 15, width = 5
pdf("ICI_projectR.TCGA.lm.SKCM.pdf", height = 12, width = 6)
p
dev.off()

# repeat this with metastatic vs non-metastatic samples
TL.proj.SKCM.nonmet <- subset(TL.proj.SKCM, definition != "Metastatic")
TL.proj.SKCM.met <- subset(TL.proj.SKCM, definition == "Metastatic")
# first for nonmets
p <- lmPlot(TL.proj.SKCM.nonmet,
            o = "OS.time",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")
# save at size height = 15, width = 5
pdf("ICI_projectR.TCGA.lm.SKCM.nonmet.pdf", height = 12, width = 6)
p
dev.off()
# now for mets
p <- lmPlot(TL.proj.SKCM.met,
            o = "OS.time",
            c = mod_cols,
            v = TRUE,
            p = "Pattern",
            f = "TL.proj$")
# save at size height = 15, width = 5
pdf("ICI_projectR.TCGA.lm.SKCM.mets.pdf", height = 12, width = 6)
p
dev.off()

# As an alternative way of visualizing mets vs nonmets
# make a plot of the delta of coefficients
lm.met <- lm(TL.proj.SKCM.met$OS.time~TL.proj.SKCM.met$Pattern_1+
              TL.proj.SKCM.met$Pattern_2+TL.proj.SKCM.met$Pattern_3+
              TL.proj.SKCM.met$Pattern_4+TL.proj.SKCM.met$Pattern_5+
              TL.proj.SKCM.met$Pattern_6+TL.proj.SKCM.met$Pattern_7+
              TL.proj.SKCM.met$Pattern_8+TL.proj.SKCM.met$Pattern_9+
              TL.proj.SKCM.met$Pattern_10+TL.proj.SKCM.met$Pattern_11+
              TL.proj.SKCM.met$Pattern_12+TL.proj.SKCM.met$Pattern_13+
              TL.proj.SKCM.met$Pattern_14+TL.proj.SKCM.met$Pattern_15+
              TL.proj.SKCM.met$Pattern_16+TL.proj.SKCM.met$Pattern_17+
              TL.proj.SKCM.met$Pattern_18+TL.proj.SKCM.met$Pattern_19+
              TL.proj.SKCM.met$Pattern_20+TL.proj.SKCM.met$Pattern_21+
              TL.proj.SKCM.met$age)
coefs.met <- summary(lm.met)$coefficients
xvals.met <- coefs.met[grep("Pattern", rownames(coefs.met)),1]
CIlow.met <- xvals.met - (1.96*coefs.met[grep("Pattern", rownames(coefs.met)),2])
CIhigh.met <- xvals.met + (1.96*coefs.met[grep("Pattern", rownames(coefs.met)),2])
# nonmet
lm.nonmet <- lm(TL.proj.SKCM.nonmet$OS.time~TL.proj.SKCM.nonmet$Pattern_1+
               TL.proj.SKCM.nonmet$Pattern_2+TL.proj.SKCM.nonmet$Pattern_3+
               TL.proj.SKCM.nonmet$Pattern_4+TL.proj.SKCM.nonmet$Pattern_5+
               TL.proj.SKCM.nonmet$Pattern_6+TL.proj.SKCM.nonmet$Pattern_7+
               TL.proj.SKCM.nonmet$Pattern_8+TL.proj.SKCM.nonmet$Pattern_9+
               TL.proj.SKCM.nonmet$Pattern_10+TL.proj.SKCM.nonmet$Pattern_11+
               TL.proj.SKCM.nonmet$Pattern_12+TL.proj.SKCM.nonmet$Pattern_13+
               TL.proj.SKCM.nonmet$Pattern_14+TL.proj.SKCM.nonmet$Pattern_15+
               TL.proj.SKCM.nonmet$Pattern_16+TL.proj.SKCM.nonmet$Pattern_17+
               TL.proj.SKCM.nonmet$Pattern_18+TL.proj.SKCM.nonmet$Pattern_19+
               TL.proj.SKCM.nonmet$Pattern_20+TL.proj.SKCM.nonmet$Pattern_21+
               TL.proj.SKCM.nonmet$age)

coefs.nonmet <- summary(lm.nonmet)$coefficients
xvals.nonmet <- coefs.nonmet[grep("Pattern", rownames(coefs.nonmet)),1]
CIlow.nonmet <- xvals.nonmet - (1.96*coefs.nonmet[grep("Pattern", rownames(coefs.nonmet)),2])
CIhigh.nonmet <- xvals.nonmet + (1.96*coefs.nonmet[grep("Pattern", rownames(coefs.nonmet)),2])
xvals.delta <- xvals.met - xvals.nonmet
# combine confidence intervals
range.delta <- sqrt((CIhigh.met - CIlow.met)**2 + (CIhigh.nonmet - CIlow.nonmet)**2)
CIlow.delta <- xvals.delta - range.delta/2
CIhigh.delta <- xvals.delta + range.delta/2
plt.df.metdelta <- data.frame(xvals = xvals.delta, 
                              CIlow = CIlow.delta,
                              CIhigh = CIhigh.delta,
                              yvals = gsub("TL.proj.SKCM.met$", "",
                                  names(xvals.delta), fixed = T))

# plot

# p <- ggplot(plt.df.metdelta,
  #       aes(y = xvals)) + 
  # geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  # geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height =
  #                  .2, color = "gray50") +
  # geom_point(col="black") +
  # geom_bar() +
  # theme_bw() +
  # #theme(panel.grid.minor = element_blank()) +
  # ylab("") +
  # xlab("Delta Standardized Coefficients")

# save at size height = 15, width = 5
pdf("ICI_projectR.TCGA.lm.SKCM.deltamets.pdf", height = 15, width = 5)
barplot(plt.df.metdelta$xvals)
dev.off()

### B7 ###

# Correlation betwen B7 and patterns in all cancers and SKCM specifically
expTCGA.t <- t(expTCGA)
expTCGA.t.B7 <- expTCGA.t[,which(colnames(expTCGA.t) == "CD80" | colnames(expTCGA.t) == "CD86")]
mtch <- match(rownames(TL.proj), rownames(expTCGA.t.B7))
TL.proj.B7 <- cbind(TL.proj, expTCGA.t.B7[mtch,])

# all TCGA samples and B7.1
patterns.B7.1.all <- apply(TL.proj.B7[,c(1:21)], 2,
              function(a){
                res <- cor.test(a, TL.proj.B7[,33])
                return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.1.df <- as.data.frame(t(patterns.B7.1.all))
colnames(patterns.B7.1.df) <- c("Estimate","CIlow",
                                 "CIhigh", "pval")
patterns.B7.1.df$pattern <- factor(rownames(patterns.B7.1.df), levels=rownames(patterns.B7.1.df))
patterns.B7.1.df$plotcol <- "black"
patterns.B7.1.df$plotcol[patterns.B7.1.df$pval <= 0.05] <- "red"
# write data to outfile
write.csv(patterns.B7.1.df, file = "ICI_projectR.B7.1.csv",
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

pdf("ICI_projectR.B7.1.Pearson.pdf", height = 9, width = 6)
p
dev.off()

# all TCGA samples and B7.2

patterns.B7.2.all <- apply(TL.proj.B7[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7[,34])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.2.df <- as.data.frame(t(patterns.B7.2.all))
colnames(patterns.B7.2.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
# if pval is zero, set to min for graphing
patterns.B7.2.df$pval[patterns.B7.2.df$pval == 0] <- min(patterns.B7.2.df$pval[patterns.B7.2.df$pval != 0])

patterns.B7.2.df$pattern <- factor(rownames(patterns.B7.2.df), levels=rownames(patterns.B7.2.df))
patterns.B7.2.df$plotcol <- "black"
patterns.B7.2.df$plotcol[patterns.B7.2.df$pval <= 0.05] <- "red"
#View(patterns.B7.2.df)
# write data to outfile
write.csv(patterns.B7.2.df, file = "ICI_projectR.B7.2.csv",
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

pdf("ICI_projectR.B7.2.Pearson.pdf", height = 9, width = 6)
p
dev.off()

# now repeat but just with SKCM
TL.proj.B7.SKCM <- TL.proj.B7[TL.proj.B7$cancertype == "TCGA-SKCM",]

# SKCM TCGA samples and B7.1
patterns.B7.1.SKCM <- apply(TL.proj.B7.SKCM[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7.SKCM[,33])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.1.SKCM.df <- as.data.frame(t(patterns.B7.1.SKCM))
colnames(patterns.B7.1.SKCM.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
patterns.B7.1.SKCM.df$pattern <- factor(rownames(patterns.B7.1.SKCM.df), levels=rownames(patterns.B7.1.SKCM.df))
patterns.B7.1.SKCM.df$plotcol <- "black"
patterns.B7.1.SKCM.df$plotcol[patterns.B7.1.SKCM.df$pval <= 0.05] <- "red"
# write data to outfile
write.csv(patterns.B7.1.SKCM.df, file = "ICI_projectR.B7.1.SKCM.csv",
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

pdf("ICI_projectR.B7.1.SKCM.Pearson.pdf", height = 9, width = 6)
p
dev.off()

# SKCM TCGA samples and B7.2

patterns.B7.2.SKCM <- apply(TL.proj.B7.SKCM[,c(1:21)], 2,
                           function(a){
                             res <- cor.test(a, TL.proj.B7.SKCM[,34])
                             return(c(res$estimate, res$conf.int, res$p.value))})
patterns.B7.2.SKCM.df <- as.data.frame(t(patterns.B7.2.SKCM))
colnames(patterns.B7.2.SKCM.df) <- c("Estimate","CIlow",
                                "CIhigh", "pval")
# if pval is zero, set to min for graphing
patterns.B7.2.SKCM.df$pval[patterns.B7.2.SKCM.df$pval == 0] <- min(patterns.B7.2.SKCM.df$pval[patterns.B7.2.SKCM.df$pval != 0])

patterns.B7.2.SKCM.df$pattern <- factor(rownames(patterns.B7.2.SKCM.df), levels=rownames(patterns.B7.2.SKCM.df))
patterns.B7.2.SKCM.df$plotcol <- "black"
patterns.B7.2.SKCM.df$plotcol[patterns.B7.2.SKCM.df$pval <= 0.05] <- "red"
#View(patterns.B7.2.df)
# write data to outfile
write.csv(patterns.B7.2.SKCM.df, file = "ICI_projectR.B7.2.SKCM.csv",
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

pdf("ICI_projectR.B7.2.SKCM.Pearson.pdf", height = 9, width = 6)
p
dev.off()

# barplots
pdf("ICI_projectR.B7.1.all.pdf")
barplot(t(patterns.B7.1.all)[,1], las = 2,
        main = "patterns.B7.1.all", ylab = "Pearson Correlation")
dev.off()

pdf("ICI_projectR.B7.2.all.pdf")
barplot(t(patterns.B7.2.all)[,1], las = 2,
        main = "patterns.B7.2.all", ylab = "Pearson Correlation")
dev.off()

pdf("ICI_projectR.B7.1.SKCM.pdf")
barplot(t(patterns.B7.1.SKCM)[,1], las = 2,
        main = "patterns.B7.1.SKCM", ylab = "Pearson Correlation")
dev.off()

pdf("ICI_projectR.B7.2.SKCM.pdf")
barplot(t(patterns.B7.2.SKCM)[,1], las = 2,
        main = "patterns.B7.2.SKCM", ylab = "Pearson Correlation")
dev.off()

# Kaplan Myers curves
TL.proj.SKCM.KM <- subset(TL.proj.SKCM, definition == "Metastatic")
Pat <- TL.proj.SKCM.KM$Pattern_7
OS <- TL.proj.SKCM.KM$OS.time
group <- rep(NA, length(OS))
event <- rep(1, length(OS))
TL.proj.km <- data.frame(OS = OS, event = event, group = group)
cutoffs <- quantile(Pat, probs = c(0.95, 0.95))
# stratify based on pattern 7
TL.proj.km$group[Pat <= cutoffs[1]] <- "Bottom95%"
TL.proj.km$group[Pat >= cutoffs[2]] <- "Top5%"
TL.proj.km <- TL.proj.km[!is.na(TL.proj.km$group),]
TL.proj.km$group <- as.factor(TL.proj.km$group)
# calc survival model fit
KMfit <- survfit(Surv(OS, event) ~ group,
                data = TL.proj.km)
# plot with survminer
pdf("ICI_projectR.TCGA.SKCM.KMplot.pdf")
ggsurvplot(KMfit, data = TL.proj.km, risk.table = TRUE, pval = T,
           palette = c("blue", "red"))
dev.off()
