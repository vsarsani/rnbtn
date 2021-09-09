
#' rnbtn_mean_agg constructs mean and control mean for the
#' nested/unnested factors listed in the fctrel from the dataframe

#' @param df : dataframe containing counts,covariates in the long format
#' @param tncnt : column corresponding to counts(y) .Ex: 'gene'
#' @param locus_tag : column corresponding to gene names/locus tags .Ex: 'gene'
#' @param fctrel : A list of column names and desired factor relevels .The order of fctrel is important.
#' First element of each covariate is taken as CONTROL




#' @import plyr
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import stringr

#' @examples
#' #Simulating and selecting Counts
#' TC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[1]]
#' #Selecting only first five hundred locus tags as an  example
#' locuslist <- TC_df$locus_tag[1:500]
#' TC_500_df <- subset(TC_df, locus_tag %in% locuslist)
#' #Preparing covariate desired levels for fct_relevel
#' fct_rel <- list(strain=c("strain_1","strain_2","strain_3"),
#' condition=c("condition_1","condition_2","condition_3","condition_4"),
#' slevel=c("slevel_1","slevel_2","slevel_3"))
#' #Calculating Control and covariate means
#' df_mean <- rnbtn_mean_agg(TC_500_df,tncnt = 'tncnt',fctrel = fct_rel)
#' @export

rnbtn_mean_agg <- function(df, tncnt = tncnt,
                           locus_tag = locus_tag, fctrel = NONE) {

  ## validate input df
  if (is.data.frame(df) == FALSE) {
    stop("Input data frame is not in the data frame
         data type. Please convert  and re-run it")
  }
  if (ncol(df) < 3) {
    stop("The data frame should have atleast one covariate,
         one response and locus/gene name column")
  }
  if (!is.null(fctrel)) {
    if (is.list(fctrel) == FALSE & is.null(names(fctrel)) == TRUE) {
      stop("Factor relevel parameter should be a list with
           names corresponding to the desired columns and levels ordered.")
    }
  }

  `%>%` <- dplyr::`%>%`

  # control mean Function
  cm <- function(df, ctrl) {
    # Take first element of every covariate list as a control
    input_list <- lapply(ctrl, `[[`, 1)
    # map with the respective covariate name
    x <- purrr::map(input_list, ~paste0(" == ", "'", .[[1]], "'"))
    # build a string of a condition for filter to parse
    conditions <- paste(names(x), x, collapse = " & ")
    # Apply the filter
    (df %>%
        dplyr::filter(eval(parse(text = conditions))))$mean
  }


  # For the list of covariates, calculate mean and control mean according to
  # structure specified in fctrel
  df_total <- data.frame()

  # for loop to capture both nested and covariate means
  for (i in 0:(length(fctrel) - 1)) {
    cond <- 1:(length(fctrel) - i)
    cond_names <- names(fctrel[cond])
    # prepare list of columns for ddply
    cond_l <- c(names(fctrel[cond]), "locus_tag")
    # ctrl input for cm function
    ctrl <- fctrel[cond_names]
    df_m <- plyr::ddply(df, cond_l, summarize, mean = mean(tncnt)) %>%
      mutate(controlmean = cm(., ctrl))
    # Combining all means while retaining original data frame
    df_total <- plyr::rbind.fill(df_total, df_m)
  }



  # Select covariates mentioned in fct_rel from the data frame
  cov_sel <- df_total %>% dplyr::select(names(fctrel))
  # Create an effect column to match model dataframe
  effect_df <- as.data.frame(sapply(colnames(cov_sel), function(name) {
    paste(name, cov_sel[, name], sep = "")
  })) %>%
    tidyr::unite(effect, 1:ncol(.), sep = ":") 
  df_total["effect"] <- effect_df
  # Prepare strings to replace 
  find.list <- paste(":",names(fct_rel), "NA", sep="")
  find.string <- paste(unlist(find.list), collapse = "|")
  df_total <- df_total%>%dplyr::mutate(effect=gsub(find.string, replacement = "", x = effect))
  ## validate output df
  if (is.data.frame(df_total) == FALSE) {
    stop("Output data frame is not in the data
         frame data type. Please check and re-run it")
  }

  return(df_total)
}
