#' Function to sample from proposed distribution
#' 
#' @export
#' @param ard matrix of respondent data
#' @param grouping vector of group for each respondent
#' @param num.iterations number of iterations to run
#' @param num.chains number of chains to run
#' @param pars parameters to save. only extract necessary rho_j by default
#' @param ... additional arguments to pass to stan like seed, cores, messages, etc.
#' @return posterior samples of rho
#' 
runSampling <- function(ard, grouping, num.iterations, 
                         num.chains, pars=c('rho_j'), ...){
    num.respondents <- nrow(ard)
    num.populations <- ncol(ard)
    num.municipalities <- length(unique(grouping))

    data <- list(
            N = num.respondents,
            K = num.populations,
            y = ard,
            jj = grouping,
            J = num.municipalities
            )

    stan.model <- stanmodels$neg_binomial_partial_pooling
    fit <- sampling(stan.model, data = data, iter = num.iterations, 
                    chains = num.chains, pars=pars, ...)
    return(fit)
}