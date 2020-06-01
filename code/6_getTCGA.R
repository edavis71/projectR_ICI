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
c <- 0
for (ct in cancer_types){
  c <- c + 1
  print(paste0("Saving and combining data from ", ct, " (", c, "/",
              length(cancer_types), ")"))
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
    tempmeta$paper_RPPA <- as.numeric(tempmeta$paper_RPPA)
  }
  if ("paper_Purity" %in% colnames(tempmeta)){
    tempmeta$paper_Purity <- as.numeric(tempmeta$paper_Purity)
  }
  if ("paper_EBV.positive" %in% colnames(tempmeta)){
    tempmeta$paper_EBV.positive <- as.numeric(tempmeta$paper_EBV.positive)
  }
  if ("paper_TP53.mutation" %in% colnames(tempmeta)){
    tempmeta$paper_TP53.mutation <- as.numeric(tempmeta$paper_TP53.mutation)
  }
  if ("paper_ABSOLUTE.Purity" %in% colnames(tempmeta)){
    tempmeta$paper_ABSOLUTE.Purity <- as.numeric(tempmeta$paper_ABSOLUTE.Purity)
  }
  if ("paper_PIK3CA.mutation" %in% colnames(tempmeta)){
    tempmeta$paper_PIK3CA.mutation <- as.numeric(tempmeta$paper_PIK3CA.mutation)
  }
  if ("paper_KRAS.mutation" %in% colnames(tempmeta)){
    tempmeta$paper_KRAS.mutation <- as.numeric(tempmeta$paper_KRAS.mutation)
  }
  if ("paper_ARID1A.mutation" %in% colnames(tempmeta)){
    tempmeta$paper_ARID1A.mutation <- as.numeric(tempmeta$paper_ARID1A.mutation)
  }
  if ("paper_purity" %in% colnames(tempmeta)){
    tempmeta$paper_purity <- as.numeric(tempmeta$paper_purity)
  }
  if ("paper_Subclonal_genome_fraction" %in% colnames(tempmeta)){
    tempmeta$paper_Subclonal_genome_fraction <- as.numeric(tempmeta$paper_Subclonal_genome_fraction)
  }
  if ("paper_age_at_diagnosis" %in% colnames(tempmeta)){
    tempmeta$paper_age_at_diagnosis <- as.numeric(tempmeta$paper_age_at_diagnosis)
  }
  if ("paper_age_at_initial_pathologic_diagnosis" %in% colnames(tempmeta)){
    tempmeta$paper_age_at_initial_pathologic_diagnosis <- as.numeric(tempmeta$paper_age_at_initial_pathologic_diagnosis)
  }
  if ("paper_ploidy" %in% colnames(tempmeta)){
    tempmeta$paper_ploidy <- as.numeric(tempmeta$paper_ploidy)
  }
  if ("paper_Ploidy" %in% colnames(tempmeta)){
    tempmeta$paper_Ploidy <- as.numeric(tempmeta$paper_Ploidy)
  }
  if ("paper_BRAF_mut" %in% colnames(tempmeta)){
    tempmeta$paper_BRAF_mut <- as.numeric(tempmeta$paper_BRAF_mut)
  }
  if ("paper_HRAS_mut" %in% colnames(tempmeta)){
    tempmeta$paper_HRAS_mut <- as.numeric(tempmeta$paper_HRAS_mut)
  }
  if ("paper_IDH1_mut" %in% colnames(tempmeta)){
    tempmeta$paper_IDH1_mut <- as.numeric(tempmeta$paper_IDH1_mut)
  }
  if ("paper_RB1_mut" %in% colnames(tempmeta)){
    tempmeta$paper_RB1_mut <- as.numeric(tempmeta$paper_RB1_mut)
  }
  if ("paper_PTEN_mut" %in% colnames(tempmeta)){
    tempmeta$paper_PTEN_mut <- as.numeric(tempmeta$paper_PTEN_mut)
  }
  if ("paper_TP53_mut" %in% colnames(tempmeta)){
    tempmeta$paper_TP53_mut <- as.numeric(tempmeta$paper_TP53_mut)
  }
  if ("paper_BAP1" %in% colnames(tempmeta)){
    tempmeta$paper_BAP1 <- as.character(tempmeta$paper_BAP1)
  }
  if ("paper_miRNA" %in% colnames(tempmeta)){
    tempmeta$paper_miRNA <- as.numeric(tempmeta$paper_miRNA)
  }
  if ("paper_Number.pack.years.smoked" %in% colnames(tempmeta)){
    tempmeta$paper_Number.pack.years.smoked <- as.numeric(tempmeta$paper_Number.pack.years.smoked)
  }
  if ("paper_Number.pack.years.smoked" %in% colnames(tempmeta)){
    tempmeta$paper_Number.pack.years.smoked <- as.numeric(tempmeta$paper_Number.pack.years.smoked)
  }
  if ("paper_Age.at.diagnosis" %in% colnames(tempmeta)){
    tempmeta$paper_Age.at.diagnosis <- as.character(tempmeta$paper_Age.at.diagnosis)
  }
  if ("paper_mRNA.cluster" %in% colnames(tempmeta)){
    tempmeta$paper_mRNA.cluster <- as.numeric(tempmeta$paper_mRNA.cluster)
  }
  # save this file
  saveRDS(tempmeta, paste0("TCGA_RDS/TCGA.legacy.meta.", ct, ".rds"))
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
  # save this file
  saveRDS(tempexp, paste0("TCGA_RDS/TCGA.legacy.expression.", ct, ".rds"))
  # add to combined meta dataframe
  if (is.null(expTCGA)){
    expTCGA <- tempexp
  } else{
    expTCGA <- cbind(expTCGA, tempexp)
  }
}

#warnings() # all should just be coercion from factors to characters

# save combined files as rds
saveRDS(metaTCGA, "TCGA_RDS/TCGA.legacy.meta.rds")
saveRDS(expTCGA, "TCGA_RDS/TCGA.legacy.expression.rds")
