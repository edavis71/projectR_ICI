#############################################################################
# Make TCGA GEP for upload to CIBERSORTx
#############################################################################

#############################################################################
# Set up environment
#############################################################################

# set working directory
setwd('/Users/michaelkessler/Dropbox/Workspace/POSTDOC/ImmunoOncology/projectR_ICI')

###################################################
### read in data from rds files I made (per cancer type if needed)
###################################################
# note these files are too big for CIBERSORTx server. Make files < 100MB each if you end up wanting to try to run all of TCGA through, one subset at a time. For now, just trying on SKCM data
# expTCGA <- readRDS("TCGA_RDS/TCGA.legacy.expression.rds")
# # split data frame (too big to run at once through CIBERSORTx)
# expTCGA1 <- expTCGA[,1:(ncol(expTCGA)/2)] 
# expTCGA2 <- expTCGA[,(ncol(expTCGA)/2):ncol(expTCGA)]
# # write outfiles
# write.table(expTCGA1, "TCGA.legacy.expression.GEP1.txt", quote = F)
# write.table(expTCGA2, "TCGA.legacy.expression.GEP2.txt", quote = F)

expTCGA <- readRDS("TCGA_RDS/TCGA.legacy.expression.TCGA-SKCM.rds")
# remove duplicate gene names
expTCGA <- expTCGA[match(unique(rownames(expTCGA)), rownames(expTCGA)),]

write.table(expTCGA, "TCGA.legacy.expression.SKCM.GEP.txt", quote = F, sep = "\t")

