#' rnbtn uses regularized negative binomial regression to estimate the
#' change in transposon insertions attributable to gene-environment
#' changes without transformations or uniform normalization
#'
#' This function constructs mean and control mean for the nested
#' factors listed in the fctrel
#'
#' @param df : dataframe containing counts,covariates in the long format
#' @param tncnt : column corresponding to counts(y) .Ex: "gene"
#' @param locus_tag : column corresponding to gene names/locus tags .Ex: "gene"
#' @param fctrel : A list of column names and desired factor relevels .
#' Example :  rlist<- list(condition = c("CONTROL","T1","T2",),
#' strain = c("WT","KO","DO"),slevel = c("none","LOW","MEDIUM","HIGH"),
#'batch = c("1","2","3","4")).
#' The order of fctrel is important.This function will assume nested order of
#' condition/strain/slevel in above example and CONTROL,WT,none as
#' respective control means for each covariate
#'
#' @examples
#'
#'tncnt="counts"
#'fctrel <- list(condition = c("CONTROL","DIS1","DIS2"),
#'strain = c("WT","KP","DO"),slevel = c("none","LOW","MEDIUM","HIGH"))
#'rnbtn_means_agg(df,tncnt,locus_tag="gene_name",fctrel=fctrel)
#'
#'
#' @export
#'
rnbtn_mean_agg <- function(df, tncnt = tncnt,
                           locus_tag = locus_tag, fctrel = NONE) {
    # Required packages tidyverse,reshape and doParallel
    suppressMessages(require(tidyverse))
    suppressMessages(require(plyr))
  ## validate input df
    if (is.data.frame(df) == FALSE) {
    stop("Input data frame is not in the data frame data type.
         Please convert  and re-run it")
  }
    if (ncol(df) < 3) {
    stop("The data frame should have atleast one covariate,
         one response and locus/gene name column")
  }
if (!is.null(fctrel)) {
    if (is.list(fctrel) == FALSE & is.null(names(fctrel)) == TRUE) {
    stop("Factor relevel parameter should be a list with names
         corresponding to the desired columns and levels ordered.")
  }
    }


# control mean Function
cm <- function(df, ctrl) {
    # Take first element of every covariate list as a control
    input_list <- lapply(ctrl, `[[`, 1)
    # map with the respective covariate name
    x <- map(input_list, ~ paste0(" == ", "'", .[[1]], "'"))
    # build a string of a condition for filter to parse
    conditions <- paste(names(x), x, collapse = " & ")
    # Apply the filter
    (df %>% filter(eval(parse(text = conditions))))$mean
    }


# For the list of covariates, calculate mean and
# control mean according to structure specified in fctrel
df_total <- data.frame()

# for loop to capture both nested and covariate means
for (i in 0 :(length(fctrel) -1)) {
    cond <- 1 :(length(fctrel)-i)
    cond_names <- names(fctrel[cond])
# prepare list of columns for ddply
cond_l <- c(names(fctrel[cond]), "locus_tag")
#ctrl input for cm function
ctrl <- fctrel[cond_names]
df_m <- ddply(df, cond_l, summarize, mean = mean(tncnt)) %>%
  mutate(controlmean = cm(., ctrl))
# Combining all means while retaining original data frame
df_total <- rbind.fill(df_total, df_m)
}



  # Select covariates mentioned in fct_rel from the data frame
cov_sel <- df_total %>% dplyr::select(names(fctrel))
# Create an effect column to match model dataframe
effect_df <- as.data.frame(sapply(colnames(cov_sel), function(name) {
  paste(name, cov_sel[, name], sep = "")})) %>%
  unite(effect, 1:ncol(.), sep = ":") %>%
  mutate_at("effect", str_replace, ":.*NA", "")
df_total["effect"] <- effect_df

 ## validate output df
    if (is.data.frame(df_total) == FALSE) {
    stop("Output data frame is not in the data frame
         data type. Please check and re-run it")
  }

return(df_total)
}
