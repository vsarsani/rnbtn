#' rnbtn uses regularized negative binomial regression to estimate the change in transposon insertions attributable to gene-environment changes 
#without transformations or uniform normalization
#' 
#' This function constructs locfdr value for each mean effect
#' 
#' @param df : dataframe of model results from  rnbtn_model_agg.R or rnbtn_model_agg_parallel.R
#' @param dfmean : dataframe of means from rnbtn_mean_agg.R
#' '@param locus_tag : column corresponding to gene names/locus tags .Ex: "gene"
#' @param exeffects : A list of effects to exclude from localfdr package
#' @param fecutoff : A cutoff for seperating frankly essential genes. Default is one. This will filter out genes/locus tags with < user defined cutoff
#' @param cecutoff: A cutoff for conditionally essential gens to be excluded from fdr calculations.Default is one. This will apply filter: pmax(mean, controlmean) > cutoff and pmin(mean, controlmean) > cutoff
#' @param path: A path where you want frankly essential and conditionally essential files written: ex: path <- "../../reports/fdr_csv/" in windows
#'@param method: Method used for fdr calculations. locfdr or fdrtool. default is locfdr
#'
#'

#' 
#' @examples 
#' 
#exeffects <- c('Intercept','batch','log(N)')
#path <- "../../reports/fdr_csv/"
#rnbtn_fdr_effects(df,dfmean,exeffects=exeffects,fecutoff=1,cecutoff=1,path,method=locfdr)

#' 
#' 
#' @export
#' 

rnbtn_fdr_effects <- function(df,dfmean,locus_tag="locus_tag",exeffects=c('Intercept','batch','log'),fecutoff=1,cecutoff=1,path,method="locfdr"){ 
    
    # Required packages tidyverse,reshape and doParallel
    suppressMessages(require(tidyverse))
    suppressMessages(require(locfdr))
    suppressMessages(require(fdrtool))
    suppressMessages(require(reshape))
    
  ## validate input df
    if(is.data.frame(df)==FALSE){
    stop("Input Model data frame is not in the data frame data type. Please convert  and re-run it")
  }
    if(ncol(df)< 3){
    stop("The data frame should have locus_tag,effect and coefficient")
  }
if(is.data.frame(dfmean)==FALSE){
    stop("Input Means data frame is not in the data frame data type. Please convert  and re-run it")
  }

    if(ncol(df)< 4){
    stop("The Means data frame should have locus_tag,effect,mean and controlmean")
  }


# make a list of all unique effects from model dataframe
effects <- as.list(unique(df$effect))
locuslist <- unique(df$locus_tag)

# Exclude effects provided by user input 
for (i in exeffects){
    effects <- effects[!grepl(i, effects)]}

# Create two empty dataframes, one for frankly essential and another one for conditionally essential
    
fe_df <- list()
ce_df <- list()
    
for ( i in effects)  {
 tryCatch({
     # apply filter for frankly essential genes
      fe_df[[i]] <- subset(dfmean, effect==i)%>%group_by(locus_tag)%>%filter(mean<fecutoff & controlmean<cecutoff)
    #extract  locus_tags and log2values
locus <- subset(df, effect == i)$locus_tag
vl <- subset(df, effect == i )$coeff_log2value
     #filter out any NA coeff values
coef.df <- data.frame(locus_tag=locus, effect=i,value=vl)%>%filter(value!="NA")
colnames(coef.df)<- c("locus_tag","effect","log2_coefficient")
count.df <- subset(dfmean, effect==i)%>%group_by(locus_tag)
coef.df <- left_join(coef.df, count.df, by=c('locus_tag', 'effect'))
#Normalize coefficient values for performing local fdr
coef.df <- coef.df%>%filter(log2_coefficient!=0) %>%mutate(n_log2_coefficient = (log2_coefficient - mean(log2_coefficient))/sd(log2_coefficient))
     
# Apply locfdr or fdrtool depending on user input
if(method=="locfdr"){
     
coef.df <- coef.df%>%filter(pmax(mean, controlmean) > cecutoff)%>%filter(pmin(mean, controlmean) > cecutoff)%>%mutate(fdr=(locfdr(n_log2_coefficient,plot=0)$fdr))%>%dplyr::select(-n_log2_coefficient)
    
     }
else
    {
 coef.df <- coef.df%>%filter(pmax(mean, controlmean) > cecutoff)%>%filter(pmin(mean, controlmean) > cecutoff)%>%mutate(fdr=(fdrtool(n_log2_coefficient,cutoff.method="locfdr",plot=FALSE,verbose=FALSE)$lfdr))%>%dplyr::select(-n_log2_coefficient) 
}

  if (is.null(coef.df)) {    
      
    cat("Unable to fdr for effect ",i,"Please change fdr method or input df","\n")
      
      }
     
# gather genes/locus_tags which missed fdr calculations or cutoffs
missedresults<- list()
suppressWarnings(suppressMessages(for (g in locuslist) {
     if(nrow(coef.df%>%filter(locus_tag==g))==0 & nrow(fe_df[[i]]%>%filter(locus_tag==g))==0 ){
        
       missed <- melt(unique(coef.df$effect))
       missedresults[[g]] <- missed%>%mutate(locus_tag=g,effect=value,log2_coefficient="NA")%>%dplyr::select(-value)   
    }   
}))
# Transforming dataframe and removing dummy variable
suppressWarnings(suppressMessages(missed.df <- melt(missedresults)%>%dplyr::select(-L1)))
missed.df <- left_join(missed.df, count.df, by=c('locus_tag', 'effect'))%>%mutate(fdr="NA")
# Joining missed locus tags with others
cond_essential<- rbind(coef.df,missed.df)
# sort by fdr
ce_df[[i]]<- arrange(cond_essential,fdr)
cat(" fdr for effect ",i,"is completed","\n")
    }, error=function(e){})
    }


# Writing outputs to paths 

    
#frankly essential
suppressMessages(lapply(1:length(fe_df), function(i) write.csv(fe_df[[i]], file.path(path,paste0(str_replace_all(names(fe_df[i]),":","_"), "_essential.csv")),
        
                                    row.names = FALSE)))
# conditionally essential
                        
suppressMessages(lapply(1:length(ce_df), function(i) write.csv(ce_df[[i]], file.path(path,paste0(str_replace_all(names(ce_df[i]),":","_"), "_effect.csv")),
                                    row.names = FALSE)))
}
