###################################################
# PLEASE NOTE: this is an example run of CoGAPS
# subset for NK cells in control and aCTLA-4 treatment 
# groups for high variance genes to reduce run time
###################################################

###################################################
### init
###################################################
library(CoGAPS)
library(monocle3)
load("monocle3_myeloid_cds_updated.rda")

###################################################
### pseudotime
###################################################

nk <- cds[,pData(cds)$assigned_cell_type == "NK"]
nk <- nk[,pData(nk)$treatment == "control" | pData(nk)$treatment == "aCTLA-4"]
nk <- nk[fData(nk)$use_for_ordering == "TRUE",]
nk_matrix <- as.matrix(exprs(nk))
dim(nk_matrix)

# create new parameters object
params <- new("CogapsParams")
# view all parameters
params
params <- setParam(params, "nPatterns", 2)

result <- CoGAPS(nk_matrix, params, nIterations=50000, messages=FALSE,
                 singleCell = TRUE, sparseOptimization = TRUE, seed = 1000)


'''
parameters below were ran for the full dataset: 
-- Standard Parameters --
  nPatterns            35 
nIterations          50000 
seed                 1000 
singleCell           TRUE 
sparseOptimization   TRUE 
distributed          genome-wide 

-- Sparsity Parameters --
  alpha          0.01 
maxGibbsMass   100 

-- Distributed CoGAPS Parameters -- 
  nSets          12 
cut            35 
minNS          6 
maxNS          18 
'''
