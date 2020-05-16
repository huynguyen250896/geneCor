# geneCor v0.1.0
#### I. Introduction
The package geneSA is built to serve as a support tool for the paper "...". Automatically individually compute the Pearson's correlation coefficients of genes shared between CNA data and EXP data, and those shared between MET data and EXP data; visualize the distribution of expression of CNVcor genes and expression of METcor genes on a page; and examine the significance of each of those skewed distributions. </br> 


#### IV. Implementation
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/geneCor")
```
Call the library;
```sh
library(geneCor)
```
running example:
```sh
geneCor(cna = df1, exp1 = df2, alternative1="less", met = df3, exp2 = df4, alternative2="greater")
```
