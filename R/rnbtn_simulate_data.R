#' rnbtn_simulate_data simulates the transposon count data for Total counts
#'  and Unique Counts in a nested experimental structure of strain/condition/slevel with replicates.
#' It simulates 4000 locus_tags for each replicate.


#' @param n_strain= Number of strains or genetic backgrounds  to simulate
#' @param n_condition= Number of experimental/environmental conditions to simulate
#' for each strain
#' @param n_slevel= Number of stress levels to simulate for each condition
#' @param n_rep = Number of replicates for each stress level

#' @examples
#' rnbtn_simulate_data(n_strain=2,n_condition=2,n_slevel=2,n_rep=2)

#'@import MASS
#'@export

rnbtn_simulate_data <- function(n_strain, n_condition, n_slevel, n_rep) {

   require(MASS)

    # Simulate shape and scale parameters for the strain.hyper-parameters of
    # the Gamma distribution were drawn from uniform distributions

    strain_shape <- runif(n_strain, min = 0, max = 5)
    strain_scale <- runif(n_strain, min = 0, max = 5)

    # strain list data frames
    strain_data_TC <- list()
    strain_data_UC <- list()

    for (i in seq_len(n_strain)) {

        # shape and scale for each strain
        a <- strain_shape[i]
        b <- strain_scale[i]
        # First the dispersion parameter was sampled from a Gamma distribution
        # for each condition and for 8 intervals (l) each containing 500 genes.
        conditions <- list(rep(list(rgamma(8, a, b)), n_condition))

        # Empty lists
        condition_data_TC <- list()
        condition_data_UC <- list()


        for (j in seq_len(n_condition)) {

            slevel_data_TC <- list()
            slevel_data_UC <- list()

            for (k in seq_len(n_slevel)) {

                replicate_data_TC <- list()
                replicate_data_UC <- list()

                for (l in seq_len(n_rep)) {

                  # The number of unique insertions for each gene was sampled
                  # from a negative binomial distribution with mean parameters
                  # shared across groups of 500 genes,mu = (0:5; 1; 2; 4; 8;
                  # 16; 32; 64)

                  genes_1 <- MASS::rnegbin(500, mu = 0.5,
                                     theta = unlist(conditions[[1]][j])[1])
                  genes_2 <- MASS::rnegbin(500, mu = 1,
                                     theta = unlist(conditions[[1]][j])[2])
                  genes_3 <- MASS::rnegbin(500, mu = 2,
                                     theta = unlist(conditions[[1]][j])[3])
                  genes_4 <- MASS::rnegbin(500, mu = 4,
                                     theta = unlist(conditions[[1]][j])[4])
                  genes_5 <- MASS::rnegbin(500, mu = 8,
                                     theta = unlist(conditions[[1]][j])[5])
                  genes_6 <- MASS::rnegbin(500, mu = 16,
                                     theta = unlist(conditions[[1]][j])[6])
                  genes_7 <- MASS::rnegbin(500, mu = 32,
                                     theta = unlist(conditions[[1]][j])[7])
                  genes_8 <- MASS::rnegbin(500, mu = 64,
                                     theta = unlist(conditions[[1]][j])[8])
                  replicate_data_UC[[l]] <- c(genes_1, genes_2, genes_3,
                                              genes_4, genes_5, genes_6,
                                              genes_7, genes_8)
                  # the total transposon insertion counts were obtained by
                  # sampling from a 303 negative binomial distribution with
                  # mean = 100 and dispersion = 1 for each unique 304 insertion
                  # site previously generated
                  replicate_data_TC[[l]] <- c(genes_1, genes_2, genes_3,
                                              genes_4, genes_5, genes_6,
                                              genes_7, genes_8) * rnegbin(4000,
                                              mu = 100, theta = 1)
                }
                slevel_data_UC[[k]] <- replicate_data_UC
                slevel_data_TC[[k]] <- replicate_data_TC
            }

            condition_data_TC[[j]] <- slevel_data_TC
            condition_data_UC[[j]] <- slevel_data_UC
        }

        strain_data_TC[[i]] <- condition_data_TC
        strain_data_UC[[i]] <- condition_data_UC
    }
    # Unnest the simulated data
    TC_counts <- unlist(strain_data_TC)
    UC_counts <- unlist(strain_data_UC)

    # add column names
    Rep <- rep(rep(rep(1:n_rep, times = n_slevel,
                       each = 4000), n_condition), n_strain)
    slevel <- rep(rep(1:n_slevel, times = n_condition,
                      each = 4000 * n_rep), n_strain)
    condition <- rep(1:n_condition, times = n_slevel,
                     each = 4000 * n_rep * n_strain)
    strain <- rep(1:n_strain, each = 4000 * n_slevel * n_condition * n_rep)
    cov <- data.frame(strain, slevel, condition, Rep)
    cov_sel <- sapply(colnames(cov), function(name) {
        paste(name, cov[, name], sep = "_")
    })
    locus_tag <- rep(seq.int(4000), n_rep * n_slevel * n_condition * n_strain)
    locus_tag <- paste0("locus_tag_", locus_tag)
    # Final simulated data

    TC_data <- data.frame(locus_tag, cov_sel, tncnt = TC_counts)
    UC_data <- data.frame(locus_tag, cov_sel, tncnt = UC_counts)
    return(list(TC_data, UC_data))
}
