% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnbtn_mean_agg.R
\name{rnbtn_mean_agg}
\alias{rnbtn_mean_agg}
\title{rnbtn_mean_agg constructs mean and control mean for the
nested/unnested factors listed in the fctrel from the dataframe}
\usage{
rnbtn_mean_agg(df, tncnt = tncnt, locus_tag = locus_tag, fctrel = NONE)
}
\arguments{
\item{df}{: dataframe containing counts,covariates in the long format}

\item{tncnt}{: column corresponding to counts(y) .Ex: 'gene'}

\item{locus_tag}{: column corresponding to gene names/locus tags .Ex: 'gene'}

\item{fctrel}{: A list of column names and desired factor relevels .The order of fctrel is important.
First element of each covariate is taken as CONTROL}
}
\description{
rnbtn_mean_agg constructs mean and control mean for the
nested/unnested factors listed in the fctrel from the dataframe
}
\examples{
#Simulating and selecting Counts
TC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[1]]
#Selecting only first five hundred locus tags as an  example
locuslist <- TC_df$locus_tag[1:500]
TC_500_df <- subset(TC_df, locus_tag \%in\% locuslist)
#Preparing covariate desired levels for fct_relevel
fct_rel <- list(strain=c("strain_1","strain_2","strain_3"),
condition=c("condition_1","condition_2","condition_3","condition_4"),
slevel=c("slevel_1","slevel_2","slevel_3"))
#Calculating Control and covariate means
df_mean <- rnbtn_mean_agg(TC_500_df,tncnt = 'tncnt',fctrel = fct_rel)
}
