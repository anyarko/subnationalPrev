#' Function to sample from proposed distribution
#' 
#' @export
#' @param ard matrix of respondent data
#' @param pop.sub.regions vector of population sizes for subregions
#' @param grouping vector of group for each respondent
#' @param num.iterations number of iterations to run
#' @param num.chains number of chains to run
#' @param num.cores number of cores to use
#' @param seed random seed
#' @return posterior samples of rho and L
#' 
run.sampling <- function(ard, pop.sub.regions, grouping, num.iterations, num.cores, num.chains, seed){
    num.respondents <- nrow(ard)
    num.populations <- ncol(ard)
    num.subregions <- length(unique(grouping))

    data <- list(
            N = num.respondents,
            K = num.populations,
            y = ard,
            offset = pop.sub.regions,
            jj = grouping,
            J = num.subregions
            )

    stan.model <- stanmodels$corr_neg_binomial_partial_pooling
    fit <- sampling(stan.model, data = data, iter = num.iterations, 
                    chains = num.chains, cores = num.cores, seed=seed, pars=c('rho_j', 'L'), 
                    verbose=F, show_messages=F)
    return(fit)
}