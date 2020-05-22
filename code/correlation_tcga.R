###################################################
# Corrleation between CTLA4 expression and CIBERSORT 
# scores from TCGA data
# Authors: Ludmila Danilova and Emily Davis-Marcisak
###################################################

###################################################
### init
###################################################
library(gplots)
library(pheatmap)
library(reshape2)

###################################################
### functions
###################################################

## finds common samples using the first 15 characters from barcode
getCommonSamples = function(samples1, samples2, nameStart = 6, nameEnd = 15, uniq = T)
{
	# finds common samples
	comSamples = intersect(substr(samples1,nameStart,nameEnd),substr(samples2,nameStart,nameEnd))
#	print(comSamples)
	# finds ids of common samples for meth and expression datasets
	methSamp = vector()
	exprSamp = vector()
	for (samp in comSamples)
	{
		ind1 = grep(samp, samples1)
		ind2 = grep(samp, samples2)
		if (length(ind1)>0 && length(ind2)>0)
		{
			if (uniq)
			{
				methSamp = c(methSamp, samples1[ind1[1]])
				exprSamp = c(exprSamp, samples2[ind2[1]])
			}else{
				methSamp = c(methSamp, samples1[ind1])
				exprSamp = c(exprSamp, samples2[ind2])			
			}
		}
	}
	return (list("samp1" = methSamp,"samp2" = exprSamp))
}

###################################################
### load data 
###################################################

## read in CTLA-4 expression data
ctla4Expr = read.csv('TCGA_CTLA4expr_allCancers.csv', row.names = 1)
table(ctla4Expr$type)
types = levels(ctla4Expr$type)

## load CIBERSORT data 
## from: https://www.cell.com/immunity/fulltext/S1074-7613(18)30121-3 PMID:29628290 PMCID: PMC5982584
ciberRes = read.table(file = 'https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a', sep = '\t',header = T, stringsAsFactors = F)
## fix sample names
ciberRes[,1] = gsub('.','-',ciberRes[,1],fixed = T)
ciberRes[,1] = substr(ciberRes[,1],1,15)
## remove duplications
dup = which(!duplicated(ciberRes[,1]))
ciberRes = ciberRes[dup,]
rownames(ciberRes) = ciberRes[,1]
ciberRes = ciberRes[,-1]

###################################################
### correlation analysis
###################################################

## correlation between CTLA4 expression and CIBERSORT cell types for every cancer type
ctla4Cor_ciber = c()
for(i in types)
{
	s = ctla4Expr[(which(ctla4Expr$type == i)),]
	comm = intersect(rownames(s),rownames(ciberRes))
	ctla4Cor_ciber = cbind(ctla4Cor_ciber, cor(cbind(ctla4Expr[comm,1],ciberRes[comm,2:23]), method = 'spearman')[,1])
}
ctla4Cor_ciber = ctla4Cor_ciber[-1,]
colnames(ctla4Cor_ciber) = types

## subset for immunogenic cancers
cancer_list <- c("SKCM", "LUSC", "LUAD", "BLCA", "KIRC","KIRP")
ctla4Cor_ciber <- ctla4Cor_ciber[,colnames(ctla4Cor_ciber) %in% cancer_list]

## waterfall style bar plot for all cancers, NK cell types only.
nk_cor <- ctla4Cor_ciber[rownames(ctla4Cor_ciber) == c("NK.cells.resting", "NK.cells.activated"),]
library(reshape2)
mnk_cor <- melt(nk_cor)

ggplot(mnk_cor, aes(x=reorder(Var2, -value), y=value, fill=Var1)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Spearman correlation")


## for resting NK cells
for(i in types)
{
  s = ctla4Expr[(which(ctla4Expr$type == i)),]
  comm = intersect(rownames(s),rownames(ciberRes))
  res <- cor.test(ctla4Expr[comm,1],ciberRes[comm,12], method = "pearson")
  print(i)
  print(paste("pvalue:", res$p.value))
  print(paste("corr:", res$estimate))
}

## for active NK cells
for(i in types)
{
  s = ctla4Expr[(which(ctla4Expr$type == i)),]
  comm = intersect(rownames(s),rownames(ciberRes))
  res <- cor.test(ctla4Expr[comm,1],ciberRes[comm,13], method = "pearson")
  print(i)
  print(paste("pvalue:", res$p.value))
  print(paste("corr:", res$estimate))
}

#########################################
## boxplot of cibersort estimation
########################################

cs <- read.table(file = 'TCGA.Kallisto.fullIDs.cibersort.relative.tsv', sep = '\t', header = TRUE)
cancer_list <- c("KIRC", "SKCM", "BLCA", "LUAD", "LUSC", "KIRP")
cs_sub <- cs[cs$CancerType %in% cancer_list,]
## subsetting for NKs
cs_sub <- cs_sub[c(1,2,13,14)]
cs_melt <- melt(cs_sub)

## order factor levels
cs_melt$CancerType <- factor(cs_melt$CancerType, levels =  c("SKCM", "BLCA", "LUSC", "LUAD", "KIRP", "KIRC"))

droplevels(cs_melt$variable)
ggplot(cs_melt, aes(x=CancerType, y=value, fill = variable)) +
  geom_boxplot(outlier.shape = NA) + theme_bw() +
  theme(axis.text.x=element_text(size=rel(0.75), angle=45, hjust = 1),
        axis.ticks.x=element_blank()) + 
  scale_fill_manual(values=c("blue", "red")) + facet_wrap(~CancerType, scales = "free_x") +
  theme(panel.spacing.x=unit(0.1, "lines"),panel.spacing.y=unit(0, "lines")) +
 geom_point(aes(fill = variable, alpha = 0.8), size = 1, shape = 21, position = position_jitterdodge()) + 
  ylab("CIBERSORT scores")

