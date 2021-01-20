# -------------------------------o0o------------------------------- #
#1. PRE-PROCESSING
# -------------------------------o0o------------------------------- #

#### Library
library(CancerSubtypes)
library(tidyr)
library(tidyverse)
library(dplyr)

#### Import data
met_UM = read.table('data_methylation_hm450.txt',sep='\t',row.names = 1, check.names = FALSE, header = TRUE)
exp_UM = read.table('data_RNA_Seq_v2_mRNA_median_Zscores.txt',sep="\t",row.names = 1,  check.names = FALSE, header=TRUE)
cna_UM = read.table('data_CNA.txt', sep = '\t', row.names = 1, check.names = FALSE, header = TRUE)

####Remove Entrez_Gene_Id
exp_UM = exp_UM[,-1]
cna_UM = cna_UM[,-1]
met_UM = met_UM[,-1]

#### Retain patients sharing among three datasets
scna_UM<-colnames(cna_UM)
sexp_UM<-colnames(exp_UM)
smet_UM<-colnames(met_UM)
s_cna_exp_UM<-intersect(scna_UM,sexp_UM)
s_cna_exp_met_UM<-intersect(smet_UM,s_cna_exp_UM)

# Only retain shared samples among 3 datasets
c_cna_UM = cna_UM[,s_cna_exp_met_UM]
c_met_UM = met_UM[,s_cna_exp_met_UM]
c_exp_UM = exp_UM[,s_cna_exp_met_UM]

#### Transform list to matrix
c_exp_UM = as.matrix(c_exp_UM)
c_cna_UM = as.matrix(c_cna_UM)
c_met_UM = as.matrix(c_met_UM)

####check missing value existing?
table(is.finite(c_exp_UM))
table(is.finite(c_cna_UM)) 
table(is.finite(c_met_UM))

#### Impute missing value
# Remove probes >50% missing
c_exp_UM <- c_exp_UM[rowSums(is.na(c_exp_UM)) < (ncol(c_exp_UM) * .5), ]
c_met_UM <- c_met_UM[rowSums(is.na(c_met_UM)) < (ncol(c_met_UM) * .5), ]
# Check missing value again
table(is.finite(c_exp_UM))
table(is.finite(c_met_UM))


#### Impute using knn method provided in the function data.imputation
c_exp_UM = data.imputation(c_exp_UM, fun = "microarray")
c_met_UM = data.imputation(c_met_UM, fun = "microarray")

#### Shared probes are reserved between EXP + CNA and EXP + MET
GENE_cna_exp_UM<-intersect(rownames(c_cna_UM),rownames(c_exp_UM))
length(GENE_cna_exp_UM) 
c_cna_UM = c_cna_UM[GENE_cna_exp_UM,]
c_exp1_UM = c_exp_UM[GENE_cna_exp_UM,]
GENE_met_exp_UM<-intersect(rownames(c_met_UM),rownames(c_exp_UM))
length(GENE_met_exp_UM) 
c_met_UM = c_met_UM[GENE_met_exp_UM,]
c_exp2_UM = c_exp_UM[GENE_met_exp_UM,]

# -------------------------------o0o------------------------------- #
#2. SKEWNESS using the tool geneCor
# -------------------------------o0o------------------------------- #

#### Library

#### Correlation between corresponding EXP + CNA
x=c_cna_UM %>% t()
y=c_exp1_UM %>% t()

#### Correlation between corresponding EXP + MET
x1=c_met_UM %>% t()
y1=c_exp2_UM %>% t()

#### Tool for Identification of CNAexp and METexp, Visualization of the distribution of expression of CNAexp genes and expression of METexp genes on a page, and Examination of the significance of each of those skewed distributions.
geneCor(dat1 = x, cordat1 = y, alternative1="less", dat2 = x1, cordat2 = y1, alternative2="greater") #compute Spearman's Rank correlation coefficients.
