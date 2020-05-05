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
# TL.proj$ESTIMATE.immune.score <- as.numeric(metaTCGA$paper_ESTIMATE.immune.score)
# # add lymphocytes
# TL.proj$paper_LYMPHOCYTE.SCORE <- as.numeric(metaTCGA$paper_LYMPHOCYTE.SCORE)
# TL.proj$paper_Percent.Lymphocyte.Infiltration <- as.numeric(metaTCGA$paper_Percent.Lymphocyte.Infiltration)
# # add other lymphoctye score
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$paper_Percent.Lymphocyte.Infiltration[i])){
#     if (!is.na(metaTCGA$paper_Percent.lymphocyte.infiltration[i])){
#       TL.proj$paper_Percent.Lymphocyte.Infiltration[i] <-
#         metaTCGA$paper_Percent.lymphocyte.infiltration[i]
#     }
#   }
# }
# # add leukocytes
# TL.proj$paper_Estimated.Leukocyte.Percentage <- as.numeric(metaTCGA$paper_Estimated.Leukocyte.Percentage)
# # add other leukocyte score
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$paper_Estimated.Leukocyte.Percentage[i])){
#     if (!is.na(metaTCGA$paper_Estimated.leukocyte.percentage[i])){
#       TL.proj$paper_Estimated.Leukocyte.Percentage[i] <-
#         metaTCGA$paper_Estimated.leukocyte.percentage[i]
#     }
#   }
# }
# TL.proj$paper_Survival <- as.numeric(metaTCGA$paper_Survival)
# TL.proj$paper_Survival_Months <- as.numeric(metaTCGA$paper_Survival..months.)
# # combine days to death data
# TL.proj$days_to_death <- metaTCGA$days_to_death
# # add paper_Days.to.death
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$days_to_death[i])){
#     if (!is.na(metaTCGA$paper_Days.to.death[i])){
#       TL.proj$days_to_death[i] <-
#         metaTCGA$paper_Days.to.death[i]
#     }
#   }
# }
# # add paper_days_to_death
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$days_to_death[i])){
#     if (!is.na(metaTCGA$paper_days_to_death[i])){
#       TL.proj$days_to_death[i] <-
#         metaTCGA$paper_days_to_death[i]
#     }
#   }
# }
# # add paper_Days.to.Death
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$days_to_death[i])){
#     if (!is.na(metaTCGA$paper_Days.to.Death[i])){
#       TL.proj$days_to_death[i] <-
#         metaTCGA$paper_Days.to.Death[i]
#     }
#   }
# }
# # add paper_CLIN.days_to_death
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$days_to_death[i])){
#     if (!is.na(metaTCGA$paper_CLIN.days_to_death[i])){
#       TL.proj$days_to_death[i] <-
#         metaTCGA$paper_CLIN.days_to_death[i]
#     }
#   }
# }
# # add paper_Days.until.death
# for (i in 1:nrow(TL.proj)){
#   if (is.na(TL.proj$days_to_death[i])){
#     if (!is.na(metaTCGA$paper_Days.until.death[i])){
#       TL.proj$days_to_death[i] <-
#         metaTCGA$paper_Days.until.death[i]
#     }
#   }
# }
# # convert days_to_death to numeric
# TL.proj$days_to_death <- as.numeric(TL.proj$days_to_death)
# add age

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

# patterns x cancer type, each independently

# for each patterns and for each cancer, see if there are any signficant linear
# relationships betweens a) pattern x age b) pattern by survival measures c)
# pattern by race/gender
# cts.ptrns <- c()
# pvals <- c()
# R2s <- c()
# for (p in 1:21){
#  for (ct in cancer_types){
#    TL.proj.m.sub <- subset(TL.proj.m, variable == paste0("Pattern_", p) &
#                              cancertype == ct)
#    TL.proj.m.sub$value <- as.numeric(TL.proj.m.sub$value)
#    mod <- lm(TL.proj.m.sub$OS.time~TL.proj.m.sub$value)
#    pval <- summary(mod)$coefficients[2,4]
#    R2 <- summary(mod)$r.squared
#    # store values in lists
#    cts.ptrns <- append(cts.ptrns, paste0("Pattern_", p, " Cancer Type: ", ct))
#    pvals <- append(pvals, pval)
#    R2s <- append(R2s, R2)
#  }
# }

#mod.df <- data.frame(CancerxPattern = cts.ptrns, pval = pvals, R2 = R2s)

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

# pdf("test.lm1.pdf", height = 20, width = 20)
# avPlots(lm1)
# dev.off()
# summary(mod.overall)

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
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = "black") +
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

# plot(x_vals,
#      y_vals,
#      yaxt = "n",
#      ann = FALSE,
#      frame.plot = FALSE,
#      pch = 19,
#      cex = 1
#      )

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

# pdf("test.lm1.pdf", height = 20, width = 20)
# avPlots(lm1)
# dev.off()
# summary(mod.overall)

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
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = "black") +
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

# pdf("test.lm1.pdf", height = 20, width = 20)
# avPlots(lm1)
# dev.off()
# summary(mod.overall)

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
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = "black") +
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

pdf("ICI_projectR.TCGA.lm.age_vs_patterns.pdf")
p
dev.off()

# pattern 7 associated - PRAD, BRCA
# explore patterns within SKCM
TL.proj.SKCM <- subset(TL.proj, cancertype == "TCGA-SKCM")
lm.SKCM <- lm(TL.proj.SKCM$OS.time~TL.proj.SKCM$Pattern_7+TL.proj.SKCM$Pattern_1+
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
yvals <- gsub("TL.proj$", "", names(xvals), fixed = T)
# fix yvals order
yvals_p7 <- yvals[1]
yvals_order <- yvals[2:length(yvals)]
yvals_order <- c(yvals_order[1:6], yvals_p7, yvals_order[7:length(yvals_order)])
dat.df <- data.frame(yvals = factor(yvals, levels = yvals_order),
                     xvals = xvals,
                     CIlow = CIlow,
                     CIhigh = CIhigh,
                     pval = pval)
# make plot using ggplot
p <- ggplot(dat.df, aes(x = xvals, y = yvals, size = -log10(pval))) + 
  geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = CIhigh, xmin = CIlow), size = .5, height = 
                   .2, color = "gray50") +
  geom_point(color = "black") +
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

pdf("ICI_projectR.TCGA.lm.SKCM.pdf")
p
dev.off()

# Train random forest models for the prediction of all survival from all cancer types, patterns, and age
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

# Correlation betwen B7 and patterns in SKCM
expTCGA.t <- t(expTCGA)
expTCGA.t.B7 <- expTCGA.t[,which(colnames(expTCGA.t) == "CD80" | colnames(expTCGA.t) == "CD86")]
mtch <- match(rownames(TL.proj), rownames(expTCGA.t.B7))
TL.proj.B7 <- cbind(TL.proj, expTCGA.t.B7[mtch,])
patterns.B7.1.all <- cor(TL.proj.B7[,c(1:21),], TL.proj.B7[,32])
patterns.B7.2.all <- cor(TL.proj.B7[,c(1:21),], TL.proj.B7[,33])
TL.proj.B7.SKCM <- TL.proj.B7[TL.proj.B7$cancertype == "TCGA-SKCM",]
patterns.B7.1.SKCM <- cor(TL.proj.B7.SKCM[,c(1:21),], TL.proj.B7.SKCM[,32])
patterns.B7.2.SKCM <- cor(TL.proj.B7.SKCM[,c(1:21),], TL.proj.B7.SKCM[,33])

pdf("patterns.B7.1.all.pdf")
barplot(patterns.B7.1.all[,1], las = 2,
        main = "patterns.B7.1.all", ylab = "Pearson Correlation")
dev.off()

pdf("patterns.B7.2.all.pdf")
barplot(patterns.B7.2.all[,1], las = 2,
        main = "patterns.B7.2.all", ylab = "Pearson Correlation")
dev.off()

pdf("patterns.B7.1.SKCM.pdf")
barplot(patterns.B7.1.SKCM[,1], las = 2,
        main = "patterns.B7.1.SKCM", ylab = "Pearson Correlation")
dev.off()

pdf("patterns.B7.2.SKCM.pdf")
barplot(patterns.B7.2.SKCM[,1], las = 2,
        main = "patterns.B7.2.SKCM", ylab = "Pearson Correlation")
dev.off()


# compare xCell proportions and patterns 
# read TCGA data
TCGA.xCell <- read.table("xCell_TCGA_RSEM.txt", sep = "\t", header = T)
# save sample names
TCGA.xCell.colnames <- TCGA.xCell$X
# remove col X before normalization (contains sample names)
TCGA.xCell <- TCGA.xCell[,2:ncol(TCGA.xCell)]
# normalize values per sample
TCGA.xCell <- apply(TCGA.xCell, 2, function(x) as.numeric(as.character(x))/sum(as.numeric(as.character(x)),na.rm = T))
# traspose xCell dfs
TCGA.xCell <- t(TCGA.xCell)
TCGA.xCell <- as.data.frame(TCGA.xCell)
View(TCGA.xCell)
# only keep TCGA samples
TCGA.xCell <- TCGA.xCell[grepl("TCGA", rownames(TCGA.xCell)),]
# reassign sample names as rownames
rownames(TCGA.xCell) <- gsub(".", "-",
                             substring(rownames(TCGA.xCell), 1, 12),
                             fixed = T)
# set col names
colnames(TCGA.xCell) <- TCGA.xCell.colnames
# match tcga and xcell rows
mtch <- match(TL.proj$barcode, rownames(TCGA.xCell))
TL.proj$NK <- TCGA.xCell[mtch,]$`NK cells`
TL.proj$NKT <- TCGA.xCell[mtch,]$NKT
TL.proj$Tregs <- TCGA.xCell[mtch,]$Tregs
TL.proj$Th1 <- TCGA.xCell[mtch,]$`Th1 cells`
TL.proj.dis <- subset(TL.proj, cancertype == "TCGA-KIRC")
ggplot(TL.proj,
        aes(x = Pattern_7, y = NKT)) +
        geom_point() +
        geom_smooth(method = "lm", col = "red") +
        theme_classic() +
        facet_wrap(~cancertype, scales = "free")

# compare CIBERSORTx proportions and patterns
# get cibersort results for TCGA from supp file from Thorsson et al. 2018
TCGAthorsson <- as.data.frame(read_excel("1-s2.0-S1074761318301213-mmc2.xlsx"))
# subset and transpose
TCGAcbr <- TCGAthorsson[,37:58]
# assign TCGA ID as col name
rownames(TCGAcbr) <- TCGAthorsson$`TCGA Participant Barcode`
# recode colnames
colnames(TCGAcbr) <- sapply(sapply(colnames(TCGAcbr), function(x) strsplit(x, "...", fixed = T)), "[", 1)
View(rownames(TCGAcbr))
# match tcga and cibersort rows
mtch <- match(TL.proj$barcode, rownames(TCGAcbr))
TL.proj$NKa <- as.numeric(TCGAcbr[mtch,]$`NK Cells Activated`)
TL.proj$NKr <- as.numeric(TCGAcbr[mtch,]$`NK Cells Resting`)
TL.proj$Tregs <- as.numeric(TCGAcbr[mtch,]$`T Cells Regulatory Tregs`)
TL.proj$T8 <- as.numeric(TCGAcbr[mtch,]$`T Cells CD8`)
TL.proj$MacM2 <- as.numeric(TCGAcbr[mtch,]$`Macrophages M2`)
TL.proj.dis <- subset(TL.proj, cancertype == "TCGA-KIRC")
ggplot(TL.proj,
       aes(x = Pattern_7, y = NKr)) +
  geom_point() +
  geom_smooth(method = "lm", col = "red") +
  theme_classic() +
  facet_wrap(~cancertype, scales = "free")

  
# explore kaplan meyer survival
library(survminer)
library(RTCGA.clinical)
library(survival)

Pat <- subset(TL.proj, cancertype == "TCGA-SKCM")$Pattern_7
OS <- subset(TL.proj, cancertype == "TCGA-SKCM")$OS.time
group <- rep(NA, length(OS))
event <- rep(1, length(OS))
TL.proj.km <- data.frame(OS = OS, event = event, group = group)
cutoffs <- quantile(Pat, probs = c(0.8, 0.8))
# stratify based on pattern 7
TL.proj.km$group[Pat <= cutoffs[1]] <- "low"
TL.proj.km$group[Pat >= cutoffs[2]] <- "high"
#TL.proj.km <- TL.proj.km[!is.na(TL.proj.km$group),]
#TL.proj.km$group <- as.factor(TL.proj.km$group)
# calc survival model fit
KMfit <- survfit(Surv(OS, event) ~ group,
               data = TL.proj.km)
# plot with survminer
ggsurvplot(KMfit, data = TL.proj.km, risk.table = TRUE, pval = T)

# explore clustering and heatmap visualizations

# see if cancers cluster by patterns
col_vector = color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
anno_cols <- sample(col_vector, length(levels(TL.proj$cancertype)))
pdf("col_test.pdf")
plot(1:length(anno_cols), rep(10,length(anno_cols)), pch = 19, cex = 2,
     col = anno_cols)
dev.off()
anno_list <- list()
for (j in 1:length(levels(TL.proj$cancertype))){
  anno_list[[levels(TL.proj$cancertype)[j]]] <- anno_cols[j]
}
row_ha = rowAnnotation("Cancer Types" = TL.proj$cancertype,
                       col = list("Cancer Types" = unlist(anno_list)))
                         
                       #   anno_block(gp = gpar(fill = rainbow(length(levels(TL.proj$cancertype)))),
                       #              labels = levels(TL.proj$cancertype), 
                       #              labels_gp = gpar(col = "white",
                       #              fontsize = 10)),
                       # annotation_name_rot = 90)
# scale data
TL.proj.rownames <- rownames(TL.proj)
TL.proj.scl <- apply(TL.proj[,1:21], 2, scale)
rownames(TL.proj.scl) <- TL.proj.rownames # reassign rownames
# open connection to plotting device
tiff(file = paste("TL.all_CoGAPs_patterns.ComplexHeatmap.TCGA.tiff"),
     width = 6000, height = 5000, units = "px", res = 800)
Heatmap(TL.proj.scl, col=inferno(100), name = "mat",
        clustering_distance_rows = "pearson",
        column_names_gp = gpar(fontsize = 5),
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 1,
        show_row_names = FALSE
)
dev.off()

# now see if patterns correlate with time to death
TL.proj.dtd <- subset(TL.proj, !is.na(paper_Survival))
# set cols
anno_list <- list()
TL.proj.dtd$cancertype <- TL.proj.dtd$cancertype
for (j in 1:length(levels(TL.proj.dtd$cancertype))){
  anno_list[[levels(TL.proj.dtd$cancertype)[j]]] <- anno_cols[j]
}

row_ha = rowAnnotation("Days to Death" =
                         anno_lines(TL.proj.dtd$days_to_death,
                                    axis_param= list(direction = "reverse",
                                                     labels_rot = 90),
                                    add_points = F),
                       "Cancer Types" = TL.proj.dtd$cancertype,
                       col = list("Cancer Types" = unlist(anno_list)))

                       # "Immune Score" =
                       #   anno_lines(TL.proj$ESTIMATE.immune.score,
                       #              axis_param= list(direction = "reverse",
                       #                               labels_rot = 90),
                       #              add_points = F),

# scale data
TL.proj.dtd.rownames <- rownames(TL.proj.dtd)
TL.proj.dtd.scl <- apply(TL.proj.dtd[,1:21], 2, scale)
rownames(TL.proj.dtd.scl) <- TL.proj.dtd.rownames # reassign rownames
# open connection to plotting device
tiff(file = paste("TL.all_CoGAPs_patterns.ComplexHeatmap.TCGA.dtd.tiff"),
     width = 6000, height = 5000, units = "px", res = 800)
Heatmap(TL.proj.dtd.scl, col=inferno(100), name = "mat",
        clustering_distance_rows = "pearson",
        column_names_gp = gpar(fontsize = 5),
        #top_annotation=ha_column, 
        left_annotation=row_ha,
        border = TRUE,
        row_km = 1,
        show_row_names = FALSE
)
dev.off()

# correlation between days to death (dtd) and patterns
View(TL.proj.scl)
TL.proj.scl <- as.data.frame(TL.proj.scl)
TL.proj.scl$immune <- TL.proj$ESTIMATE.immune.score
TL.proj.scl$days_to_death <- TL.proj$days_to_death
TL.proj.scl$age <- TL.proj$age_at_index
TL.proj.scl$name <- rownames(TL.proj.scl)
TL.proj.scl$cancertype <- TL.proj$cancertype
TL.proj.scl$cols <- anno_cols[TL.proj$cancertype]
TL.proj.scl.m <- reshape2::melt(TL.proj.scl,
      id.vars = c("immune", "days_to_death","age", "cancertype", "cols", "name"),
      na.rm = FALSE, value.name = "value")

# explore each pattern
TL.proj.scl.m$value <- as.numeric(TL.proj.scl.m$value)
TL.proj.scl.m.pat <- subset(TL.proj.scl.m, variable == "Pattern_4")
pdf("TL.TCGA.pattern_vs_outcomes.pdf")
ggplot(TL.proj.scl.m.pat,
       aes(x = value, y = age, col = "black")) +
       geom_point() +
       geom_smooth(method = "lm", se = FALSE) +
       scale_colour_manual(values = anno_cols) +
       theme_classic() +
       facet_wrap(~cancertype, scales = "free")
dev.off()

