#' @title geneCor: Identification of correlation between two pairs of three datasets, visualization of Z-score distributions of 
#' those two pairs on a page, and examination of the significance of the skewness of those distributions.
#'
#' @description  It automatically computes correlation coefficients of individual genes that share between the first data and its
#' corresponding third data, and those that share between the second data and its corresponding third data; visualizes the Z-score 
#' distributions of between the first and second data versus their corresponding third data on a page; and examines the significance 
#' of the skewness for those distributions using D'Agostino test. For further information on requirements as well as how to implement 
#' this tool, please visit my Github repository: https://github.com/huynguyen250896/geneCor.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage geneCor(dat1, cordat1, alternative1, dat2, cordat2, alternative2, methodCC, adjustedP)
#'
#' @param dat1 data.frame or matrix. The first input data includes its rows are samples and its columns are genes.
#'
#' @param cordat1 data.frame or matrix. The data includes its rows are samples and its columns are clinical features. 
#' This itself is the corresponding third data of \code{dat1}. Namely, correlation analysis will be implemented 
#' between \code{dat1} and \code{cordat1}.
#'
#' @param alternative1 a character string specifying the alternative hypothesis for Z-score distribution between \code{dat1} 
#' and \code{cordat1}. Must be one of "two.sided", "greater" or "less". You can specify just the initial letter.
#'
#' @param dat2 data.frame or matrix. The second input data includes its rows are samples and its columns are genes.
#'
#' @param cordat1 data.frame or matrix. The data includes its rows are samples and its columns are clinical features. 
#' This itself is the corresponding third data of \code{dat2}. Namely, correlation analysis will be implemented 
#' between \code{dat2} and \code{cordat2}.
#'
#' @param alternative2 a character string specifying the alternative hypothesis for Z-score distribution between \code{dat2} 
#' and \code{cordat2}. Must be one of "two.sided", "greater" or "less". You can specify just the initial letter.
#'
#' @param methodCC character. correlation method. Allowed values are \code{spearman}, \code{pearson} (default), \code{kendall}.
#'
#' @param adjustedP logical. Whether we should adjust the P-values gained from correlation analyses. Default is \code{adjustedP = T}
#'
#' @return NULL
#'
#' #compute Spearman's Rank correlation coefficients.
#' @examples geneCor(dat1 = cna, cordat1 = exp1, alternative1="less", dat2 = met, cordat2 = exp2, alternative2="greater", method = "spearman") 
#' #dat1 receives copy number alterations data, and cordat1 receives its corresponding gene expression data.
#' #dat2 receives methylation data, and cordat2 receives its corresponding gene expression data.
#'
#' @export

geneCor = function(dat1 = NULL, cordat1 = NULL, alternative1=c("two.sided","less","greater"), dat2 = NULL, cordat2 = NULL, alternative2=c("two.sided","less","greater"), methodCC = "pearson", adjustedP = TRUE)
{
  #library
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(moments)
  
  #ERRORs
  if(missing(dat1)){
    stop("Error: dat1 is missing \n")
  }
  
  if(missing(cordat1)){
    stop("Error: cordat1 is missing \n")
  }
  
  if(all(dim(dat1) != dim(cordat1))){
    stop("Error: dimension of dat1 must equal dimension of cordat1 to perform correlation analysis between them \n")
  }
  
  if(missing(dat2)){
    stop("Error: dat2 is missing \n")
  }
  
  if(missing(cordat2)){
    stop("Error: cordat2 is missing \n")
  }
  
  if(all(dim(dat2) != dim(cordat2))){
    stop("Error: dimension of dat2 must equal dimension of cordat2 to perform correlation analysis between them \n")
  }
  
  if(missing(alternative1)){
    stop("Error: Please specify the alternative hypothesis for Z-score distribution between dat1 and cordat1 \n")
  }
  
  if(missing(alternative2)){
    stop("Error: Please specify the alternative hypothesis for Z-score distribution between dat2 and cordat2 \n")
  }
  
  if(all(rownames(dat1) != rownames(cordat1))){
    stop("Error: patients in dat1 must be included in cordat1 at their rows and in exactly the same order \n")
  }

  if(all(rownames(dat2) != rownames(cordat2))){
    stop("Error: patients in dat2 must be included in cordat2 at their rows and in exactly the same order \n")
  }
  
  if(all(colnames(dat1) != colnames(cordat1))){
    stop("Error: genes in dat1 must be included in cordat1 at their columns and in exactly the same order \n")
  }
  
  if(all(colnames(dat2) != colnames(cordat2))){
    stop("Error: genes in dat2 must be included in cordat2 at their columns and in exactly the same order \n")
  }
  
  #define the computeQ function may adjust gained P-values following Benjamini-Hochberg FDR
  computeQ <- function(x)
  {
    (x$P.value*nrow(x))/(x$rank)
  }
  
  #dat1
  dat1_cor<-data.frame(My_name_is=paste("Huy", 1:ncol(dat1)),CC=NA ,P.value=NA)
  estimates = numeric(ncol(dat1))
  pvalues = numeric(ncol(dat1))
  for (i in c(1:ncol(dat1))) {
    cc=cor.test(dat1[,i],cordat1[,i], method = methodCC)
    dat1_cor$CC[i]=cc$estimate
    dat1_cor$P.value[i]=cc$p.value
    rownames(dat1_cor) = colnames(dat1)[1:ncol(dat1)]
  }
  dat1_cor = dat1_cor[,-1]
  if(adjustedP == TRUE | adjustedP == T){
    order.pvalue = order(dat1_cor$P.value)
    dat1_cor = dat1_cor[order.pvalue,] #order rows following p-value
    dat1_cor$rank = rank(dat1_cor$P.value) #re-order
    dat1_cor$Q.value = computeQ(dat1_cor) #compute Q-value
    dat1_cor = dat1_cor %>% subset(P.value <= 0.05) #only retain Genes with P <=0.05
    dat1_cor = dat1_cor %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
    dat1_cor = dplyr::select(dat1_cor, -rank) #remove the 'rank' column
  } else{
    dat1_cor= dat1_cor %>% subset(P.value <=0.05)
  }
  dat1_cor$fisher_z_trans = log((1+dat1_cor$CC)/(1-dat1_cor$CC))
  
  #dat2
  dat2_cor <- data.frame(My_name_is=paste("Huy", 1:ncol(dat2)),CC=NA ,P.value=NA)
  estimates = numeric(ncol(dat2))
  pvalues = numeric(ncol(dat2))
  for (i in c(1:ncol(dat2))) {
    cc1=cor.test(dat2[,i],cordat2[,i], method = methodCC)
    dat2_cor$CC[i]=cc1$estimate
    dat2_cor$P.value[i]=cc1$p.value
    rownames(dat2_cor) = colnames(dat2)[1:ncol(dat2)]
  }
  dat2_cor = dat2_cor[,-1]
  if(adjustedP == T | adjustedP == TRUE){
    order.pvalue = order(dat2_cor$P.value)
    dat2_cor = dat2_cor[order.pvalue,] #order rows following p-value
    dat2_cor$rank = rank(dat2_cor$P.value) #re-order
    dat2_cor$Q.value = computeQ(dat2_cor) #compute Q-value
    dat2_cor = dat2_cor %>% subset(P.value <= 0.05) #only retain Genes with P <=0.05
    dat2_cor = dat2_cor %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
    dat2_cor = dplyr::select(dat2_cor, -rank) #remove the 'rank' column
  } else{
  dat2_cor = dat2_cor %>% subset(P.value <=0.05)}
  dat2_cor$fisher_z_trans = log((1+dat2_cor$CC)/(1-dat2_cor$CC))
  cat("- Correlation analyses are performed successfully...", "\n")
  cat("- Gained correlation coefficients are converted to Z values by Fisher’s Z-transformation following the equation: Z = 0.5*ln[(1+r)/(1−r)].", "\n")
  
  #### p-value of skewness using D'Agostino skewness test
  cat("- Examine whether outliers exist in Z-score values between dat1 and cordat1 or not...", "\n")
  cat(">> The number of negative Z-scores is", table(sign(dat1_cor$fisher_z_trans))[[1]], ", whereas the number of positive Z-scores is", table(sign(dat1_cor$fisher_z_trans))[[2]],"\n")
  if(table(sign(dat1_cor$fisher_z_trans))[[1]] > table(sign(dat1_cor$fisher_z_trans))[[2]]){
    answer = as.character(readline('Perhaps positive Z-scores are potential outliers. Do you want to remove them? (Y or N?) \n'))
    if(answer == "Y"){
      dat1_cor = dat1_cor %>% dplyr::filter(.$fisher_z_trans < 0)
    }
  } else{
    answer = as.character(readline('Perhaps negative Z-scores are potential outliers. Do you want to remove them? (Y or N?) \n'))
    if(answer == "Y"){
      dat1_cor = dat1_cor %>%  dplyr::filter(.$fisher_z_trans > 0)
    }
  }
  
  cat("- Examine whether outliers exist in Z-score values between dat2 and cordat2 or not...", "\n")
  cat(">> The number of negative Z-score is", table(sign(dat2_cor$fisher_z_trans))[[1]], ", whereas the number of positive Z-score is", table(sign(dat2_cor$fisher_z_trans))[[2]],"\n")
  if(table(sign(dat2_cor$fisher_z_trans))[[1]] > table(sign(dat2_cor$fisher_z_trans))[[2]]){
    answer = as.character(readline('Perhaps positive Z-scores are potential outliers. Do you want to remove them? (Y or N?) \n'))
    if(answer == "Y"){
      dat2_cor = dat2_cor %>% dplyr::filter(.$fisher_z_trans < 0)
    }
  } else{
    answer = as.character(readline('Perhaps negative Z-scores are potential outliers. Do you want to remove them? (Y or N?) \n'))
    if(answer == "Y"){
      dat2_cor = dat2_cor %>%  dplyr::filter(.$fisher_z_trans > 0)
    }
  }
  
  set.seed(25896)
  print(agostino.test(dat1_cor$fisher_z_trans, alternative = alternative1))
  print(agostino.test(dat2_cor$fisher_z_trans, alternative = alternative2))
  
  ####  visualization of distribution
  cat("- Print Z-score distribution of between dat1 versus cordat1, and dat2 versus cordat2...", "\n")
  dat1_name = as.character(readline('Please name dat1 explicitly. What name you want to set in abbreviation for dat1 is  (e.g., gene expression ~ GE): \n'))
  dat2_name = as.character(readline('Please name dat2 explicitly. What name you want to set in abbreviation for dat2 is (e.g., gene expression ~ GE): \n'))
  
  DF1 = as.data.frame(dat1_cor) %>% mutate(Dataset = dat1_name)
  DF2 = as.data.frame(dat2_cor) %>% mutate(Dataset = dat2_name)
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
  
  #write results in the wroking directory
  write.table(dat1_cor,"dat1_versus_cordat1.txt",sep = "\t", quote = FALSE)
  write.table(dat2_cor,"dat2_versus_cordat2.txt",sep = "\t", quote = FALSE)
  
  #warning
  writeLines("NOTE:\n*dat1_versus_cordat1.txt and dat2_versus_cordat2.txt placed in your current working directory.")
}
