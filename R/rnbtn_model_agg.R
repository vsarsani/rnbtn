
#' rnbtn_model_agg aggregates the regularized nested negative binomial
#' model results for all genes in a serial fashion.

#' @param df : Dataframe containing counts, covariates in the long format
#' @param formula : Provide model matrix formula using as.formula()
#' @param locus_tag : A column corresponding to gene names/locus tags .Ex: 'gene'
#' @param fctrel : A list of column names and desired factor relevels .
#' @param a : elastic net mixing parameter . Default is zero
#' @param iter: Number of times to run cross validation to take
#' the mean error associated with each lambda value, and
#'  then choose lambda.Default is 5. iter increases your run time


#' @import MASS
#' @import glmnet
#' @import tidyverse
#' @import dplyr
#' @import reshape
#' @import forcats

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
#' #Run and aggregrate model results
#' rnbtn_model_agg(TC_20_df,formula = formula,locus_tag = "locus_tag",
#' fctrel = fct_rel, iter =2, a=0)
#' @export

rnbtn_model_agg <- function(df, formula, locus_tag = locus_tag,
                            fctrel = NONE, iter = 5, a = 0) {


    ## validate input df
    if (is.data.frame(df) == FALSE) {
        stop("Input data frame is not in the data frame data type.
             Please convert it")
    }
    if (ncol(df) < 3) {
        stop("The data frame should have atleast one covariate,
             response and gene column")
    }
    if (!is.null(fctrel)) {
        if (is.list(fctrel) == FALSE & is.null(names(fctrel)) == TRUE) {
            stop("Factor relevel parameter should be a list with names
                 corresponding to the desired columns and levels")
        }
    }

    # Extracting Response(counts) column from formula
    fctrel <- fctrel
    res_name <- gsub("[()]", "", formula[2])
    # Checking if response is counts
    if (is.numeric(df[[res_name]]) == FALSE) {
        stop("The response specified is incorrect or  is not numeric")
    }

    #### Modeling Function
    rnbtn_model_pergene <- function(df, formula, res_name,
                                    locus_tag = locus_tag,
        iter = iter, a = a) {
        # Extract the genename from the model
        locus_tag <- as.character(unique(df$locus_tag))
        x <- model.matrix(as.formula(formula), df)
        ## Model starts here
        model <- tryCatch({
            x <- model.matrix(as.formula(formula), df)
            y <- df[[res_name]]
            lambda_seq <- 10^seq(2, -2, by = -0.1)
            t <- MASS::glm.nb(as.formula(formula), df)$theta
            lambdas <- NULL
            for (i in 1:iter) {
                fit <- glmnet::cv.glmnet(x, y, alpha = a, lambda = lambda_seq,
                                 family = negative.binomial(theta = t),
                  parallel = FALSE, intercept = FALSE)
                errors <- data.frame(fit$lambda, fit$cvm)
                lambdas <- rbind(lambdas, errors)
            }
            # take mean cvm for each lambda
            lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

            # select the best one
            bestindex <- which(lambdas[2] == min(lambdas[2]))
            best_lambda <- lambdas[bestindex, 1]
            best_ridge <- glmnet::glmnet(x, y, alpha = a, lambda = best_lambda,
                                 family = negative.binomial(theta = t),
                                 intercept = FALSE)
        }, error = function(err) {
            status <<- err$message
            return(NULL)
        })

        col_length <- ncol(x)

        ## Writing row of NA's if model is null. Usually Coefficient Output has
        ## 8 columns. Change this
        if (is.null(coef(model))) {
            cat("Unable to fit model for this locus_tag :", locus_tag, "\n")
            namatrix <- t(as.matrix(rep(NA, col_length)))
            row.names(namatrix) <- locus_tag
            return(namatrix)
        } else {
            vals <- t(as.matrix(coef(model)))
            coeffs <- t(as.matrix(vals[, -(which(colSums(vals) == 0))]))
            row.names(coeffs) <- locus_tag
            cat("model complete for this locus_tag :", locus_tag, "\n")
            return(coeffs)
        }
    }


    `%>%` <- dplyr::`%>%`
    #### Gene/locus wise Fits

    # For Each gene run the model and store results
    locusresults <- list()
    locuslist <- unique(df[[locus_tag]])
    cat("Running model in serial.This might take a while.
        If you want to speed it up, use parallel option",
        "\n")

    suppressWarnings(suppressMessages(for (i in 1:length(locuslist)) {
        g <- locuslist[i]
        # Selecting gene/locus tag from list
        df_g <- df %>% dplyr::filter(locus_tag == g)
        # Applying factor relevels
        if (!is.null(fctrel)) {
            for (i in names(fctrel)) {
                df_g[[i]] <- forcats::fct_relevel(factor(df_g[[i]]), unname(fctrel[i]))
            }
        }

        # using rnbtn_model_pergene function

        locusresults[[g]] <- suppressWarnings(suppressMessages(
          rnbtn_model_pergene(df_g,
            formula, res_name, locus_tag = "g", iter = 5, a = 0)))
        # garbage collector
        gc()

    }))
    ## Post-Process and storing results in a data frame Here L1 is dummy header
    ## given by melt function
    suppressWarnings(mod_data <- reshape::melt(locusresults) %>%
        dplyr::select(-L1))
    # Select only gene,effect and coefficient
    colnames(mod_data) <- c("locus_tag", "effect", "coeff")
    col_length <- length(unique(mod_data[["effect"]]))
    # Stop if none of models fit
    if (all(is.na(mod_data$coeff))== TRUE) {
        stop(" glmnet failed on all locus_tags. Please recheck
your dataframe or formula ")
    }



    # Scale coefficient to log2
    suppressWarnings(mod_data <- mod_data %>%
        dplyr::mutate(coeff_log2value = as.numeric(coeff) / log(2)) %>%
        dplyr::filter(!effect %in% (1:col_length)))

    # Restore genes/locus_tags that failed model as 'NA' so that we have an
    # record of them
    missedresults <- list()
    suppressWarnings(suppressMessages(for (g in locuslist) {

        if (nrow(mod_data %>%
            filter(locus_tag == g)) == 0) {

            missed <- reshape::melt(unique(mod_data$effect))
            missedresults[[g]] <- missed %>%
                dplyr::mutate(locus_tag = g, effect = value, coeff = "NA",
                       coeff_log2value = "NA") %>%
                dplyr::select(-value)
        }
    }))
    # Here L1 is dummy header given by melt function.Remove that header
    suppressWarnings(suppressMessages(missed_data <- reshape::melt(missedresults) %>%
        dplyr::select(-L1)))
    # Combining both missed and model results
    total_model_data <- rbind(mod_data, missed_data)


    # Check if your final output gene numbers are same as the inpu
    initial_tags <- length(locuslist)
    final_tags <- length(unique(total_model_data$locus_tag))

    if (initial_tags != final_tags) {
        cat("The gene/locus tags in the output are
            not equal to the initial input. ")
    } else {
        cat("The gene/locus tags in the output are
            equal to the initial input.Process Finished ")
    }

    return(total_model_data)
}
