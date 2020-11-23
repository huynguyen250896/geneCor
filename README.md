# geneCor v0.1.1
#### I. Introduction
---
The package geneCor is built to serve as a support tool for the paper "*[Multi-omics analysis detects novel prognostic subgroups of breast cancer](https://www.frontiersin.org/articles/10.3389/fgene.2020.574661/full?utm_source=F-NTF&utm_medium=EMLX&utm_campaign=PRD_FEOPS_20170000_ARTICLE#F5)*". </br> It automatically individually computes correlation coefficients of genes shared between CNA data and the corresponding mRNA, and those shared between MET data and the corresponding mRNA; visualizes the overall distribution of Z values between MET or CNA and the corresponding mRNA on a page; and examines the significance of the skewness for those distributions using D'Agostino test. </br> 

#### II. Data Struture
---
You must preprare the four kinds of the following data: *df1*, *df2*, *df3*, and *df4* (see the 'III.Implementation' section).</br>  
df1: copy-number alteration matrix whose rows are samples, and columns are genes. </br>  
df2: gene expression matrix is corresponding to the matrix df1, whose rows are samples, and columns are genes. NOTE that the size of the matrix df2 is the same as that of the matrix df1. </br>  
df3: methylation matrix whose rows are samples, and columns are genes. </br>  
df4: gene expression matrix is corresponding to the matrix df3, whose rows are samples, and columns are genes. NOTE that the size of the matrix df4 is the same as that of the matrix df3. </br>  
Please download datasets [Dataset](https://github.com/huynguyen250896/geneCor/tree/master/Dataset) as examples to well grasp geneCor's requirement on data structure. </br>

#### III. Implementation
---
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
geneCor(cna = df1, exp1 = df2, alternative1="less", met = df3, exp2 = df4, alternative2="greater") #compute Pearson's correlation coefficients (default method).
geneCor(cna = df1, exp1 = df2, alternative1="less", met = df3, exp2 = df4, alternative2="greater", method = "spearman") #compute Spearman's Rank correlation coefficients.
geneCor(cna = df1, exp1 = df2, alternative1="less", met = df3, exp2 = df4, alternative2="greater", method = "kendall") #compute Kendall's correlation coefficients.
```
#### IV. What's new
---
- 2020-09-30 : The function now can compute one of the three common correlation methods: Pearson, Spearman's rank, or Kendall's tau-b.

#### V. Citation 
---
Please kindly cite the following paper (and Star this Github repository if you find this tool of interest) if you use the tool in this repo: </br>
```sh
Author: Nguyen, Quang-Huy
Nguyen, Hung
Nguyen, Tin
Le, Duc-Hau
Year: 2020
Title: Multi-omics analysis detects novel prognostic subgroups of breast cancer
Journal: Frontiers in Genetics
Type of Article: ORIGINAL RESEARCH
DOI: 10.3389/fgene.2020.574661
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
