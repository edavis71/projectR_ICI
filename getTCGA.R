#############################################################################
# Use the TCGAbiolinks R package to download all TCGA expression data
# Store this data as R objects for future use
# Store per cancer type and combined
# Authors: Michael D. Kessler and Emily Davis Marcisak
#############################################################################

#############################################################################
# Set up environment
#############################################################################

# load packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
# set working directory
setwd('/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI')

###################################################
### download data per cancer type
###################################################
# get valid project IDs
project_IDs <- TCGAbiolinks:::getGDCprojects()$project_id
# get "TCGA" named projects
cancer_types <- project_IDs[grep("TCGA", TCGAbiolinks:::getGDCprojects()$project_id)]
# add a few other interesting projects
cancer_types <- c(cancer_types)
# download and parse data per cancer type using TCGABioLinks
TCGA_data <- list() # list of tissue specific TCGA data frames
for (ct in cancer_types){
  # use TCGABioLinks to download data and parse into dataframes
  # define data you want to download
  query <- GDCquery(project = ct,
                    data.category = "Gene expression",
                    data.type = "Gene expression quantification",
                    platform = "Illumina HiSeq", 
                    file.type  = "normalized_results",
                    experimental.strategy = "RNA-Seq",
                    legacy = TRUE)
  # download data
  GDCdownload(query, method = "api", files.per.chunk = 10)
  # parse/prepare the data 
  TCGAdata <- GDCprepare(query)
  # store cancer type specific meta data
  metaTCGA <- as.data.frame(SummarizedExperiment::colData(TCGAdata))
  # store cancer type specific expression data
  expTCGA <- as.matrix(SummarizedExperiment::assay(TCGAdata))
  # store this data in the TCGA_data list
  TCGA_data[[ct]] <- list()
  TCGA_data[[ct]][["meta"]] <- metaTCGA
  TCGA_data[[ct]][["exp"]] <- expTCGA
}

#############################################################################
# save TCGA R data structures for future use
# combine cancer type specific meta and expression data frames
# save combined data as R data structure for future use
#############################################################################
metaTCGA <- NULL # reinit
expTCGA <- NULL # reinit
for (ct in cancer_types){
  print(ct)
  # extract cancer speicfic metadata
  tempmeta <- TCGA_data[[ct]][["meta"]]
  tempmeta$cancertype <- rep(ct, nrow(tempmeta))
  
  # fix a few cols manually that are leading to row_bind errors due to type
  if ("paper_miRNA.cluster" %in% colnames(tempmeta)){
    tempmeta$paper_miRNA.cluster <- as.numeric(tempmeta$paper_miRNA.cluster)
  }
  if ("paper_Genome.doublings" %in% colnames(tempmeta)){
    tempmeta$paper_Genome.doublings <- as.numeric(tempmeta$paper_Genome.doublings)
  }
  if ("paper_Batch" %in% colnames(tempmeta)){
    tempmeta$paper_Batch <- as.numeric(tempmeta$paper_Batch)
  }
  if ("paper_RPPA" %in% colnames(tempmeta)){
    tempmeta$paper_Batch <- as.numeric(tempmeta$paper_RPPA)
  }
  
  # save this file
  saveRDS(tempmeta, paste0("TCGA_RDS/TCGA.legacy.expression.", ct, ".rds"))
  # add to combined meta dataframe
  if (is.null(metaTCGA)){
    metaTCGA <- tempmeta
  } else{
    # first find col intersect
    #colintersect <- intersect(colnames(tempmeta), colnames(metaTCGA))
    metaTCGA <- dplyr::bind_rows(metaTCGA, tempmeta) # retain all cols
  }
  # extract cancer speicfic expression data
  tempexp <- TCGA_data[[ct]][["exp"]]
  tempexp$cancertype <- rep(ct, nrow(tempexp))
  # save this file
  saveRDS(tempexp, paste0("TCGA_RDS/TCGA.legacy.meta.", ct, ".rds"))
  # add to combined meta dataframe
  if (is.null(expTCGA)){
    expTCGA <- tempexp
  } else{
    expTCGA <- cbind(expTCGA, tempexp)
  }
}

# save combined files as rds
saveRDS(metaTCGA, "TCGA_RDS/TCGA.legacy.meta.rds")
saveRDS(expTCGA, "TCGA_RDS/TCGA.legacy.expression.rds")

View(TCGA_data[["TCGA-LGG"]][["meta"]]$paper_Batch)
View(metaTCGA$paper_Genome.doublings)

  