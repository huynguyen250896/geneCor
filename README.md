# geneCor v0.1.0
#### I. Introduction
The package geneCor is built to serve as a support tool for the paper "*Multi-omics analysis detects novel prognostic subgroups of breast cancer*". </br> It automatically individually computes the Pearson's correlation coefficients of genes shared between CNA data and the corresponding mRNA, and those shared between MET data and the corresponding mRNA; visualizes the overall distribution of Z values between MET or CNA and the corresponding mRNA on a page; and examines the significance of the skewness for those distributions using D'Agostino test. </br> 

#### II. Data Struture
You must preprare the four kinds of the following data: *df1*, *df2*, *df3*, and *df4* (see the 'III.Implementation' section).</br>  
df1: . </br>  
df2: . </br>  
df3: . </br>  
df4: . </br>  
Please download datasets [Dataset](https://github.com/huynguyen250896/geneCor/tree/master/Dataset) as examples to well grasp geneCor's requirement on data structure. </br>

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

#### V. Citation 
Please kindly cite the following paper if you use the code, datasets or any results in this repo: </br>
```sh
...
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
