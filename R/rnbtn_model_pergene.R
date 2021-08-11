#' rnbtn uses regularized negative binomial regression to estimate the change in transposon insertions attributable to gene-environment changes 
#without transformations or uniform normalization
#' 
#' This function/module handles modeling part of rnbtn per gene
#' 
#' @param x : model matrix/design matrix of the experiment  
#' @param y : counts coressponsing to the experiment in matrix format. If in vector or list frmat, convert by suing as.matrix
#' @param locus_tag : the name/id of locus_tag for modeling
#' @param intercept: TRUE/FALSE. Whether you want intercept to be included in the model. Default is FALSE
#' @param parallel: TRUE/FALSE. Whether you want your cross validation steps in run in parallel. Default is TRUE
#' @param alpha : elastic net mixing parameter . Default is zero
#' @param iter: Number of times to run cross validation to take the mean error associated with each lambda value, and then choose lambda.Default is 10
# Note: iter increases your run time 
#' 
#' @examples 
#' 
#' 
#' 
#' @export
#' 


rnbtn_model_pergene <- function(x, y, locus_tag=NULL, intercept=FALSE, parallel=FALSE, 
                       alpha=0, iter=10){ 
    
    # Required packages glmnet and MASS
    suppressMessages(require(glmnet))
    suppressMessages(require(MASS))
    
  ## validate input x and y  
    if(is.matrix(x)==FALSE){
    stop("Model matrix is not in a matrix data type. Please convert it by as.matrix and re-run it")
  }
    if(is.matrix(y)==FALSE){
    stop("y(response) is not in a matrix data type. Please convert it by as.matrix and re-run it")
  }
if(nrow(x)!=nrow(y)){
    stop("Sample size lengths do not match between x(model matrix) and y(response).Please correct it and re-run")
  }

    
# Fit the model by performing regularization 
    
  model = tryCatch(
                {
                lambda_seq <- 10^seq(2, -2, by = -.1)
                # Estimate Theta value for the cv.glmnet
                suppressWarnings(suppressMessages(t <- glm.nb(y~x)$theta))
                # create empty vector of LAMBDAS
                lambdas = NULL
                
                # For each iteration of cross-validation
                for (i in 1:iter)
                   {
                if( parallel & intercept){
                cat("Running model in parallel,iteration: ",i ,"\n")
 suppressWarnings(suppressMessages(fit <- cv.glmnet(x, y, alpha = alpha, lambda = lambda_seq,family = negative.binomial(theta = t),parallel=TRUE,intercept = TRUE)))
} else if (parallel & !intercept) {
         cat("Running model in parallel,iteration:",i ,"\n")
suppressWarnings(suppressMessages(fit <- cv.glmnet(x, y, alpha = alpha, lambda = lambda_seq,family = negative.binomial(theta = t),parallel=TRUE,intercept = FALSE)))
} else if (!parallel & intercept){
                     cat("Running model without parallel,iteration:",i,"\n")
suppressWarnings(suppressMessages(fit <- cv.glmnet(x, y, alpha = alpha, lambda = lambda_seq,family = negative.binomial(theta = t),parallel=FALSE,intercept = TRUE)))
}else {cat("Running model without parallel,iteration:",i,"\n")
                    suppressWarnings(suppressMessages(fit <- cv.glmnet(x, y, alpha = alpha, lambda = lambda_seq,family = negative.binomial(theta = t),parallel=FALSE,intercept = FALSE)))}
                 # store errors and lambdas
                 errors = data.frame(fit$lambda,fit$cvm)
                 lambdas <- rbind(lambdas,errors)
                }
              # take mean cvm for each lambda
              lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)

             # select the best one
            bestindex = which(lambdas[2]==min(lambdas[2]))
            best_lambda = lambdas[bestindex,1]
            suppressWarnings(suppressMessages(best_ridge <- glmnet(x, y, alpha = 0, lambda = best_lambda,family = negative.binomial(theta = t),intercept = FALSE)))
                },
                error=function(err) {
                  status <<- err$message
                  return(NULL)
                })
    

## Writing row of NA's if model is null. Usually Coefficient Output has roughly same columns as model matrix. 
            col_length=ncol(x)
             if (is.null(coef(model))) { 
                 cat("Unable to fit model for this locus_tag \n")
                 namatrix=t(as.matrix(rep("NA",col_length)))
                 row.names(namatrix)= locus_tag
                 return (namatrix)}
            ## Writing  coefficients 
             else {
                 vals<- t(as.matrix(coef(model)))
                 coeffs=t(as.matrix(vals[,-(which(colSums(vals) == 0))]))
                  row.names(coeffs)= locus_tag
                 return (coeffs)
                 }
            }
