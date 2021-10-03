#'rnbtn_fdr_effects constructs the locfdr values for each nested/un nested effect from the model fits.

#' @param df : Dataframe of model results from  rnbtn_model_agg.R
#'  or rnbtn_model_agg_parallel.R
#' @param dfmean : Dataframe of mean/controlmean from rnbtn_mean_agg.R
#' @param locus_tag : Column corresponding to gene names/locus tags .
#'  Ex: 'gene'
#' @param exeffects : A list of effects to exclude from locfdr calculations.
#' @param fecutoff : A cutoff for separating frankly essential genes.
#' Default is one.
#' This will filter out genes/locus tags with < user defined cutoff
#' @param cecutoff: A cutoff for conditionally essential genes
#' to be excluded from fdr calculations.
#' Default is one. This will apply filter:
#' pmax(mean, controlmean) > cutoff and
#' pmin(mean, controlmean) > cutoff
#' @param method: Method used for fdr calculations.
#' locfdr or fdrtool. Default is locfdr
#'

#' @import dplyr
#' @import reshape
#' @import locfdr
#' @import fdrtool

#' @examples
#' #Simulating and selecting Counts
#' TC_df <- rnbtn_simulate_data(n_strain=3,n_condition=4,n_slevel=3,n_rep=2)[[1]]
#' #Selecting only first twenty locus tags as an  example
#' locuslist <- TC_df$locus_tag[1:20]
#' TC_20_df <- subset(TC_df, locus_tag %in% locuslist)
#' #Preparing covariate desired levels for fct_relevel
#' fct_rel <- list(strain=c("strain_1","strain_2","strain_3"),
#' condition=c("condition_1","condition_2","condition_3","condition_4"),
#' slevel=c("slevel_1","slevel_2","slevel_3"))
#' #Model nested formula
#' formula <- as.formula(tncnt ~ strain/condition/slevel)
#' #Run and aggregrate model results in parallel fashion
#' model_df <- rnbtn_model_agg_parallel(TC_20_df,formula = formula,
#' locus_tag = "locus_tag",fctrel = fct_rel,
#' iter =2, a=0, cores=2,ctype= "PSOCK")
#' #Calculating Control and covariate means
#' df_mean <- rnbtn_mean_agg(TC_20_df,tncnt = 'tncnt',fctrel = fct_rel)
#' #Calculate local fdr
#' TC_fdr <- rnbtn_fdr_effects(model_df,dfmean = df_mean,locus_tag = "locus_tag",
#' exeffects = c("Intercept"),fecutoff = 1,cecutoff = 1,method = "fdrtool")
#' @export

rnbtn_fdr_effects <- function(df, dfmean, locus_tag = "locus_tag", exeffects = c("Intercept",
    "batch", "log"), fecutoff = 1, cecutoff = 1, method = "locfdr") {


    ## validate input df
    if (is.data.frame(df) == FALSE) {
        stop("Input Model data frame is not in the data frame data type. Please convert  and re-run it")
    }
    if (ncol(df) < 3) {
        stop("The data frame should have locus_tag,effect and coefficient")
    }
    if (is.data.frame(dfmean) == FALSE) {
        stop("Input Means data frame is not in the data frame data type. Please convert  and re-run it")
    }

    if (ncol(df) < 4) {
        stop("The Means data frame should have locus_tag,effect,mean and controlmean")
    }


    # make a list of all unique effects from model dataframe
    effects <- as.list(unique(df$effect))
    locuslist <- unique(df$locus_tag)

    # Exclude effects provided by user input
    for (i in exeffects) {
        effects <- effects[!grepl(i, effects)]
    }

    # Create two empty dataframes, one for frankly essential and another one
    # for conditionally essential

    fe_df <- list()
    ce_df <- list()

    `%>%` <- dplyr::`%>%`

    for (i in effects) {
        tryCatch({
            # apply filter for frankly essential genes
            fe_df[[i]] <- subset(dfmean, effect == i) %>%
                dplyr::group_by(locus_tag) %>%
                dplyr::filter(mean < fecutoff & controlmean < cecutoff)
            # extract locus_tags and log2values
            locus <- subset(df, effect == i)$locus_tag
            vl <- subset(df, effect == i)$coeff_log2value
            # filter out any NA coeff values
            coef.df <- data.frame(locus_tag = locus, effect = i, value = vl) %>%
                filter(value != "NA")
            colnames(coef.df) <- c("locus_tag", "effect", "log2_coefficient")
            count.df <- subset(dfmean, effect == i) %>%
                dplyr::group_by(locus_tag)
            coef.df <- dplyr::left_join(coef.df, count.df, by = c("locus_tag", "effect"))
            # Normalize coefficient values for performing local fdr
            coef.df <- coef.df %>%
                dplyr::filter(log2_coefficient != 0) %>%
                dplyr::mutate(n_log2_coefficient = (log2_coefficient - mean(log2_coefficient))/sd(log2_coefficient))

            # Apply locfdr or fdrtool depending on user input
            if (method == "locfdr") {

                coef.df <- coef.df %>%
                  dplyr::filter(pmax(mean, controlmean) > cecutoff) %>%
                  dplyr::mutate(fdr = (locfdr::locfdr(n_log2_coefficient, plot = 0)$fdr)) %>%
                  dplyr::select(-n_log2_coefficient)

            } else {
                coef.df <- coef.df %>%
                  dplyr::filter(pmax(mean, controlmean) > cecutoff) %>%
                  dplyr::mutate(fdr = (fdrtool::fdrtool(n_log2_coefficient, cutoff.method = "locfdr",
                    plot = FALSE, verbose = FALSE)$lfdr)) %>%
                  dplyr::select(-n_log2_coefficient)
            }

            if (is.null(coef.df)) {

                cat("Unable to fdr for effect ", i, "Please change fdr method or input df",
                  "\n")

            }

            # gather genes/locus_tags which missed fdr calculations or cutoffs
            missedresults <- list()
            suppressWarnings(suppressMessages(for (g in locuslist) {
                if (nrow(coef.df %>%
                  dplyr::filter(locus_tag == g)) == 0 & nrow(fe_df[[i]] %>%
                  dplyr::filter(locus_tag == g)) == 0) {

                  missed <- reshape::melt(unique(coef.df$effect))
                  missedresults[[g]] <- missed %>%
                    dplyr::mutate(locus_tag = g, effect = value, log2_coefficient = "NA") %>%
                    dplyr::select(-value)
                }
            }))
            # Transforming dataframe and removing dummy variable
            suppressWarnings(suppressMessages(missed.df <- melt(missedresults) %>%
                dplyr::select(-L1)))
            missed.df <- dplyr::left_join(missed.df, count.df, by = c("locus_tag", "effect")) %>%
                dplyr::mutate(fdr = "NA")
            # Joining missed locus tags with others
            cond_essential <- rbind(coef.df, missed.df)
            # sort by fdr
            ce_df[[i]] <- dplyr::arrange(cond_essential, fdr)
            cat(" fdr for effect ", i, "is completed", "\n")
        }, error = function(e) {
        })
    }


return(list(fe_df,ce_df))
}
