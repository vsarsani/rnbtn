### About rnbtn Package

rnbtn is a model-based method that uses regularized negative binomial regression to estimate the change in transposon insertions attributable to gene-environment changes without transformations or uniform normalization. An empirical Bayes model for estimating the local false discovery rate combines unique and total count information to test for genes that show a statistically significant change in transposon counts. When applied to RB-TnSeq (randomized barcode transposon sequencing) and Tn-seq (transposon sequencing) libraries made in strains of Caulobacter crescentus using both total and unique count data the model was able to identify a set of conditionally essential genes for each target condition that shed light on their functions and roles during various stress conditions.


Our approach for integrating all of the experimental data to estimate the effect of the many genetic backgrounds and the environmental conditions is based on a generalized linear 
model framework.Estimating the model parameters when the number of transposon count is small has been noted by others and handled either by filtration 
addition of pseudo-counts . The low counts in response variables can result in inflated regression coefficients and are susceptible to very high variance. They also affect false discovery rate procedures increasing the risk of type-I errors. Hencem we apply a regularized negative binomial regression in a nexted GLM fashion.



rnbtn is available at biorxiv:

>Vishal Sarsani, Berent Aldikacti, Shai He, Rilee Zeinert, Peter Chien, Patrick Flaherty ,
Model-based identification of conditionally-essential genes from transposon-insertion sequencing data ,
doi: https://doi.org/10.1101/2021.07.15.452443

## Installation


# Install from GitHub

To use rnbtn from my GitHub, you will need to properly install devtools package, and download the rnbtbn package with

```r
# install.packages("devtools")
devtools::install_github("vsarsani/rnbtn")
```



# Install from Bioconductor
We are under the process of submitting/reviewing this package to Bioconductor.It will be available in the future release of Bioconductor.(3.14 or higher)



# Usage
Basic usage (users need to input at least three parameters - a data frame of counts and covariates, design formula and factor relevel test.)

```r
library(rnbtn)
```

Use your own dataframe or use simulated data using our rnbtn_simulate_data function.

```r
# Simulated Total counts of Tn-seq Experiment
TC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[1]]
# Simulated Unique counts of Tn-seq Experiment
UC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[2]]
```

Prepare the experimental design formula and factor list.
Note: First element in factor list of every covariate is taken as control


```r
#Preparing covariate desired levels for fct_relevel
fct_rel <- list(strain=c("strain_1","strain_2","strain_3"),
 condition=c("condition_1","condition_2","condition_3","condition_4"),
 slevel=c("slevel_1","slevel_2","slevel_3"))
 #Model nested formula
 formula <- as.formula(tncnt ~ strain/condition/slevel)
 ```


Run the regularized negative binomial regression model on all locus_tags in serial fashion

```r
model_df <- rnbtn_model_agg(TC_df,formula = formula,locus_tag = "locus_tag", fctrel = fct_rel, iter =2, a=0)
 ```

Run the regularized negative binomial regression model on all locus_tags in parallel fashion with 10 cores

```r
model_df <- rnbtn_model_agg_parallel(TC_df,formula = formula, locus_tag = "locus_tag",fctrel = fct_rel,
 iter =2, a=0, cores=10,ctype= "PSOCK")
  ```

Construct mean and control mean for the nested/unnested factors listed in the fctrel from the dataframe

```r
df_mean <- rnbtn_mean_agg(TC_df,tncnt = 'tncnt',fctrel = fct_rel)
```


Run the local fdr for each nested/unnested effect. fecutoff is cutoff for frankly essential. cecutoff is cutoff for conditionally essential 
```r
TC_fdr <- rnbtn_fdr_effects(model_df,dfmean = df_mean,locus_tag = "locus_tag",
exeffects = c("Intercept"),fecutoff = 1,cecutoff = 1,method = "fdrtool")
```



Optional: Gather all Total counts,Unique counts effect files and identify significant locus tags through intersection and plots.

```r
rnbtn_intersect_effects(TC_files, UC_files, path, fdrcutoff = 0.05)
```

