# geneCor v0.1.1
#### I. Introduction
---
The package geneCor is built to serve as a support tool for the paper "*[Multi-omics analysis detects novel prognostic subgroups of breast cancer](https://www.frontiersin.org/articles/10.3389/fgene.2020.574661/full?utm_source=F-NTF&utm_medium=EMLX&utm_campaign=PRD_FEOPS_20170000_ARTICLE#F5)*". </br> It automatically computes correlation coefficients of individual genes that share between the first data and its corresponding third data, and those that share between the second data and its corresponding third data; visualizes the Z-score distributions of between the first and second data versus their corresponding third data on a page; and examines the significance of the skewness for those distributions using D'Agostino test.  </br> 

#### II. Understanding the tool
---
The following are parameters provided by geneCor:
- dat1: data.frame or matrix. The first input data includes its rows are samples and its columns are genes.

- cordat1: data.frame or matrix. The data includes its rows are samples and its columns are clinical features. This itself is the corresponding third data of `dat1`. Namely, correlation analysis will be implemented between `dat1` and `cordat1`.

- alternative1: a character string specifying the alternative hypothesis for Z-score distribution between `dat1` and `cordat1`. Must be one of "two.sided", "greater" or "less". You can specify just the initial letter.

- dat2: data.frame or matrix. The second input data includes its rows are samples and its columns are genes.

- cordat2: data.frame or matrix. The data includes its rows are samples and its columns are clinical features. This itself is the corresponding third data of `dat2`. Namely, correlation analysis will be implemented between `dat2` and `cordat2`.

- alternative2: a character string specifying the alternative hypothesis for Z-score distribution between `dat2` and `cordat2`. Must be one of "two.sided", "greater" or "less". You can specify just the initial letter.

- methodCC: character. correlation method. Allowed values are `spearman`, `pearson` (default), `kendall`.

- adjustedP: logical. Whether we should adjust the P-values gained from correlation analyses using the Benjamini-Hochberg procedure. Default is `adjustedP = T`.

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
geneCor(dat1 = cna, cordat1 = exp1, alternative1="less", dat2 = met, cordat2 = exp2, alternative2="greater") #compute Pearson's correlation coefficients.
#' #dat1 receives copy number alterations data, and cordat1 receives its corresponding gene expression data.
#' #dat2 receives methylation data, and cordat2 receives its corresponding gene expression data.

geneCor(dat1 = cna, cordat1 = exp1, alternative1="less", dat2 = met, cordat2 = exp2, alternative2="greater", method = "spearman")  #compute Spearman's Rank correlation coefficients.

geneCor(dat1 = cna, cordat1 = exp1, alternative1="less", dat2 = met, cordat2 = exp2, alternative2="greater", method = "kendall") #compute Kendall's correlation coefficients.
```
#### IV. What's new
---
- 2021-01-20 : the function now can adjust gained P-values from the process of correlation analysis using the Benjamini-Hochberg procedure.
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
