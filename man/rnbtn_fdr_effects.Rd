% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnbtn_fdr_effects.R
\name{rnbtn_fdr_effects}
\alias{rnbtn_fdr_effects}
\title{rnbtn_fdr_effects constructs the locfdr values for each nested/un nested effect from the model fits.}
\usage{
rnbtn_fdr_effects(
  df,
  dfmean,
  locus_tag = "locus_tag",
  exeffects = c("Intercept", "batch", "log"),
  fecutoff = 1,
  cecutoff = 1,
  method = "locfdr"
)
}
\arguments{
\item{df}{: Dataframe of model results from  rnbtn_model_agg.R
or rnbtn_model_agg_parallel.R}

\item{dfmean}{: Dataframe of mean/controlmean from rnbtn_mean_agg.R}

\item{locus_tag}{: Column corresponding to gene names/locus tags .
Ex: 'gene'}

\item{exeffects}{: A list of effects to exclude from locfdr calculations.}

\item{fecutoff}{: A cutoff for separating frankly essential genes.
Default is one.
This will filter out genes/locus tags with < user defined cutoff}

\item{cecutoff:}{A cutoff for conditionally essential genes
to be excluded from fdr calculations.
Default is one. This will apply filter:
pmax(mean, controlmean) > cutoff and
pmin(mean, controlmean) > cutoff}

\item{method:}{Method used for fdr calculations.
locfdr or fdrtool. Default is locfdr}
}
\description{
rnbtn_fdr_effects constructs the locfdr values for each nested/un nested effect from the model fits.
}
\examples{
#Simulating and selecting Counts
TC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[1]]
#Selecting only first twenty locus tags as an  example
locuslist <- TC_df$locus_tag[1:20]
TC_20_df <- subset(TC_df, locus_tag \%in\% locuslist)
#Preparing covariate desired levels for fct_relevel
fct_rel <- list(strain=c("strain_1","strain_2","strain_3"),
condition=c("condition_1","condition_2","condition_3","condition_4"),
slevel=c("slevel_1","slevel_2","slevel_3"))
#Model nested formula
formula <- as.formula(tncnt ~ strain/condition/slevel)
#Run and aggregrate model results in parallel fashion
model_df <- rnbtn_model_agg_parallel(TC_20_df,formula = formula,
locus_tag = "locus_tag",fctrel = fct_rel,
iter =2, a=0, cores=2,ctype= "PSOCK")
#Calculating Control and covariate means
df_mean <- rnbtn_mean_agg(TC_20_df,tncnt = 'tncnt',fctrel = fct_rel)
#Calculate local fdr
TC_fdr <- rnbtn_fdr_effects(model_df,dfmean = df_mean,locus_tag = "locus_tag",
exeffects = c("Intercept"),fecutoff = 1,cecutoff = 1,method = "fdrtool")
}
