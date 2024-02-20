#' Function to scale posterior samples of rho by prior known group sizes
#'
#' @export
#' @param rho matrix from posterior samples of rho_j
#' @param correlation matrix from posterior samples of L
#' @param raw_known_sizes vector of prior known group sizes for subregion
#' @param population_size total population size for subregion
#' @return matrix of scaled posterior samples of rho_j
#'
correlatedScaling <- function(rho, correlation, 
                              raw_known_sizes,
                              population_size){
  iterations <- nrow(rho)
  K <- ncol(rho)
  new_rho <- rho
  known_ind <- which(raw_known_sizes > 0)
  
  for(k in 1:K){
    for(iter in 1:iterations){
      l_known_ind <- setdiff(known_ind, k)
      num_known <- length(l_known_ind)
      w <- correlation[iter, k, l_known_ind]
      w[w < 0] <- 0 
      w <- w * ( num_known / sum(w) )
      
      known_sizes <- raw_known_sizes[-k]
      known_sizes <- known_sizes[which(known_sizes > 0)]
      
      C_m  <- log( (1/(num_known)) * 
                     sum( (exp(rho[iter, l_known_ind]) * w ) / (unlist(known_sizes)/population_size)) )
      new_rho[iter, k] <- rho[iter, k] - C_m
      
    }
  }
  return(new_rho)
}
