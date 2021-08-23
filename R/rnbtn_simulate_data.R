#' rnbtn uses regularized negative binomial regression to estimate the change in transposon insertions attributable to gene-environment changes 
#without transformations or uniform normalization
#' 
#' This function simulates  transposon count data for Total counts and Unique Counts
#' 
#' @param 

#' nstrain= Number of strains/genetic backgrounds or different organisms to simulate
#' ncondition= Number of experimental or environmental conditions to simulate
#' nslevel= Number of stress levels to simulate for each condition
#' nrep = Number of replicates for each stress level 

#' 
#' @examples 
#' 
#rnbtn_simulate_data(nstrain=2,ncondition=2,nslevel=2,nrep=2)

#' 
#' 
#' @export
#' 

rnbtn_simulate_data <- function(nstrain=3,ncondition=4,nslevel=3,nrep=3,path){ 
suppressMessages(require(MASS))
    
    

# Simulate shape and scale parameters for the strain.hyper-parameters of the Gamma distribution were drawn from uniform distributions
    
strain_shape <- runif(nstrain,min=0,max=5)
strain_scale <- runif(nstrain,min=0,max=5)
    
#strain list data frames
strain_data_TC<-list()
strain_data_UC<- list()
    
for (i in 1:nstrain){

# shape and scale for each strain
a <- strain_shape[i]
b <- strain_scale[i]
#First the dispersion parameter was sampled from a Gamma distribution for each condition and for 8 intervals (l) each containing 500 genes.
conditions <- list(rep(list(rgamma(8,a,b)),ncondition))

# Empty lists
condition_data_TC<- list()
condition_data_UC<- list()

    
for (j in 1:ncondition){

slevel_data_TC<- list()
slevel_data_UC<- list()
    
for (k in 1:nslevel){

replicate_data_TC<-list()
replicate_data_UC<-list()

 for (l in 1:nrep){

#The number of unique insertions for each gene was sampled from a negative binomial distribution 
#with mean parameters shared across groups of 500 genes,mu = (0:5; 1; 2; 4; 8; 16; 32; 64)

genes_1<- rnegbin(500, mu = 0.5,theta=unlist(conditions[[1]][j])[1])
genes_2<- rnegbin(500, mu = 1,theta=unlist(conditions[[1]][j])[2])
genes_3<- rnegbin(500, mu = 2,theta=unlist(conditions[[1]][j])[3])
genes_4<- rnegbin(500, mu = 4,theta=unlist(conditions[[1]][j])[4])
genes_5<- rnegbin(500, mu = 8,theta=unlist(conditions[[1]][j])[5])
genes_6<- rnegbin(500, mu = 16,theta=unlist(conditions[[1]][j])[6])
genes_7<- rnegbin(500, mu = 32,theta=unlist(conditions[[1]][j])[7])
genes_8<- rnegbin(500, mu = 64,theta=unlist(conditions[[1]][j])[8])
replicate_data_UC[[l]] <- c(genes_1,genes_2,genes_3,genes_4,genes_5,genes_6,genes_7,genes_8)
#the total transposon insertion counts were obtained by sampling from a 303
#negative binomial distribution with mean  = 100 and dispersion  = 1 for each unique 304
#insertion site previously generated
replicate_data_TC[[l]] <- c(genes_1,genes_2,genes_3,genes_4,genes_5,genes_6,genes_7,genes_8)*rnegbin(4000,mu=100,theta=1) 
}
slevel_data_UC[[k]] <- replicate_data_UC
slevel_data_TC[[k]] <- replicate_data_TC
}

condition_data_TC[[j]] <- slevel_data_TC
condition_data_UC[[j]] <- slevel_data_UC
}

strain_data_TC[[i]]<-condition_data_TC
strain_data_UC[[i]]<-condition_data_UC
}
# Unnest the simulated data
TC_counts <- unlist(strain_data_TC)
UC_counts <- unlist(strain_data_UC)
                      
# add column names 
Rep <- rep(rep(rep(1:nrep, times=nslevel, each=4000),ncondition),nstrain)
Slevel <- rep(rep(1:nrep, times=nslevel, each=4000*ncondition),nstrain)
Condition <- rep(1:nrep, times=nslevel, each=4000*ncondition*nstrain)
Strain <- rep(1:nrep, each=4000*nslevel*ncondition*nstrain)
cov <- data.frame(Strain,Slevel,Condition,Rep)
cov_sel <- sapply(colnames(cov),function(name){ paste(name,cov_sel[,name],sep="_")})
locus_tag <- rep(seq.int(4000),nrep*nslevel*ncondition*nstrain)
locus_tag<- paste0("locus_tag_",locus_tag)
# Final simulated data
                      
TC_data <- data.frame(locus_tag,Strain,Condition,Slevel,Rep,tncnt=TC_counts)
UC_data <- data.frame(locus_tag,Strain,Condition,Slevel,Rep,tncnt=UC_counts)
return(list(TC_data,UC_data))
}
