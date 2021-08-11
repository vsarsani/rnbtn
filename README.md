About rnbtn

rnbtn is a model-based method that uses regularized negative binomial regression to estimate the change in transposon insertions attributable to gene-environment changes without transformations or uniform normalization. An empirical Bayes model for estimating the local false discovery rate combines unique and total count information to test for genes that show a statistically significant change in transposon counts. When applied to RB-TnSeq (randomized barcode transposon sequencing) and Tn-seq (transposon sequencing) libraries made in strains of Caulobacter crescentus using both total and unique count data the model was able to identify a set of conditionally essential genes for each target condition that shed light on their functions and roles during various stress conditions.


Our approach for integrating all of the experimental data to estimate the effect of the many genetic backgrounds and the environmental conditions is based on a generalized linear 
model framework.Estimating the model parameters when the number of transposon count is small has been noted by others and handled either by filtration 
addition of pseudo-counts . The low counts in response variables can result in inflated regression coefficients and are susceptible to very high variance. They also affect false discovery rate procedures increasing the risk of type-I errors. Hencem we apply a regularized negative binomial regression in a nexted GLM fashion.



rnbtn is available at biorxiv:

Model-based identification of conditionally-essential genes from transposon-insertion sequencing data
 Vishal Sarsani, Berent Aldikacti, Shai He, Rilee Zeinert, Peter Chien, Patrick Flaherty
doi: https://doi.org/10.1101/2021.07.15.452443

Installation

Install from Bioconductor
ComBat-seq is available in Bioconductor sva v3.36.0[2,3], please download and install following instructions in the link to Bioconductor sva. I am not a maintainer for the Bioconductor sva package, and may not be able to address some issues in that version (many of the issues are related to package dependencies and versions).

Install from GitHub
To use rnbtn from my GitHub, you will need to properly install devtools package, and download the rnbtbn package with

Install from Bioconductor
We are under the process of submitting/reviewing this package to Bioconductor.It might be available in the future release of Bioconductor.

# install.packages("devtools")
devtools::install_github("vsarsani/rnbtn")


Usage
Basic usage (users need to input at least two parameters - a data frame of counts, covariates and design formula.):

