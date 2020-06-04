# geneCor v0.1.0
#### I. Introduction
The package geneCor is built to serve as a support tool for the paper "...". It automatically individually computes the Pearson's correlation coefficients of genes shared between CNA data and the corresponding mRNA, and those shared between MET data and the corresponding mRNA; visualizes the overall distribution of Z values between MET or CNA and the corresponding mRNA on a page; and examines the significance of the skewness for those distributions using D'Agostino test. </br> 


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

#### V.Citation 
Please kindly cite the two repositories if you use the code, datasets or any results in this repo: </br>
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3872667.svg)](https://doi.org/10.5281/zenodo.3872667)
```sh
@software{nguyen_quang_huy_2020_3872667,
  author       = {Nguyen, Quang-Huy},
  title        = {huynguyen250896/geneCor: GeneCor v0.1.0},
  month        = jun,
  year         = 2020,
  publisher    = {Zenodo},
  version      = {v0.1.0},
  doi          = {10.5281/zenodo.3872667},
  url          = {https://doi.org/10.5281/zenodo.3872667}
}
```
Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
