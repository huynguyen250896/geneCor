#' @title Tool for Identification of CNAcor and METcor, Visualization of the distribution of expression of CNVcor genes and expression of METcor genes on a page, and Examination of the significance of the skewness of those distributions.
#'
#' @description  It automatically individually computes the Pearson's correlation coefficients of genes shared between CNA data and the corresponding mRNA, and those shared between MET data and the corresponding mRNA; visualizes the overall distribution of Z values between MET or CNA and the corresponding mRNA on a page; and examines the significance of the skewness for those distributions using D'Agostino test.
#'
#' @param cna,exp1,alternative1,met,exp2,alternative2
#'
#' @return NULL
#'
#' @examples geneCor(cna = df1, exp1 = df2, alternative1="less", met = df3, exp2 = df4, alternative2="greater")
#'
#' @export

# identify significantly correlated genes
geneCor = function(cna = NULL, exp1 = NULL, alternative1=c("two.sided","less","greater"), met = NULL, exp2 = NULL, alternative2=c("two.sided","less","greater"), method = "pearson")
{
  #library
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(moments)

  #missing inputs
  if(missing(cna)){
    stop("Error: cna is missing \n")
  }

  if(missing(exp1)){
    stop("Error: exp1 is missing \n")
  }

  if(missing(met)){
    stop("Error: met is missing \n")
  }

  if(missing(exp2)){
    stop("Error: exp2 is missing \n")
  }

  if(missing(alternative1)){
    stop("Error: alternative1 is missing \n")
  }

  if(missing(alternative2)){
    stop("Error: alternative2 is missing \n")
  }
  #CNA
  cna_cor<-data.frame(name=paste("Site", 1:ncol(cna)),Estimate=NA ,P.value=NA)
  estimates = numeric(ncol(cna))
  pvalues = numeric(ncol(cna))
  for (i in c(1:ncol(cna))) {
    cc=cor.test(cna[,i],exp1[,i], method = method)
    cna_cor$Estimate[i]=cc$estimate
    cna_cor$P.value[i]=cc$p.value
    rownames(cna_cor) = colnames(cna)[1:ncol(cna)]
  }
  cna_cor = cna_cor[,-1]
  cna_cor= cna_cor %>% subset(P.value <=0.05)
  cna_cor$fisher_z_trans = log((1+cna_cor$Estimate)/(1-cna_cor$Estimate))
  write.table(cna_cor,"cna_cor.txt",sep = "\t", quote = FALSE)

  #MET
  met_cor <- data.frame(name=paste("Site", 1:ncol(met)),Estimate=NA ,P.value=NA)
  estimates = numeric(ncol(met))
  pvalues = numeric(ncol(met))
  for (i in c(1:ncol(met))) {
    cc2=cor.test(met[,i],exp2[,i], method = method)
    met_cor$Estimate[i]=cc2$estimate
    met_cor$P.value[i]=cc2$p.value
    rownames(met_cor) = colnames(met)[1:ncol(met)]
  }
  met_cor = met_cor[,-1]
  met_cor = met_cor %>% subset(P.value <=0.05)
  met_cor$fisher_z_trans = log((1+met_cor$Estimate)/(1-met_cor$Estimate))
  write.table(met_cor,"met_cor.txt",sep = "\t", quote = FALSE)

  #### p-value of skewness
  set.seed(25896)
  print(agostino.test(cna_cor$fisher_z_trans, alternative = alternative1))

  print(agostino.test(met_cor$fisher_z_trans, alternative = alternative2))

  ####  visualization of distrubution
  DF1 = as.data.frame(cna_cor) %>% mutate(Dataset = "CNA")
  DF2 = as.data.frame(met_cor) %>% mutate(Dataset = "MET")
  DF <- bind_rows(DF1,DF2)

  DF_mean <- DF %>% group_by(Dataset) %>%
    summarise(Mean = mean(fisher_z_trans),
              Median = median(fisher_z_trans),
              Var = var(fisher_z_trans))

  q = ggplot(DF,aes(x = fisher_z_trans, fill = Dataset))+
    geom_density(alpha = 0.6)+
    geom_vline(inherit.aes = FALSE,
               data = DF_mean, aes(xintercept = Mean, color = Dataset),
               linetype = "dashed", size = 2,
               show.legend = FALSE) + xlab("Z value") 
  print(q)

  #warning
  writeLines("NOTE:\n*cna_cor.txt and met_cor.txt placed in your current working directory\n*Please check to identify which significantly DNA copy-number-correlated genes found between CNA and EXP, and which significantly methylation-correlated genes found between MET and EXP")
}
