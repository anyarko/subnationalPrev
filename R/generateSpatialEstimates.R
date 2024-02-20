#' Function to calculate the average correlation between groups
#' 
#' @export
#' @param params posterior samples
#' @return a matrix of the correlation between groups at each iteration
#' 
compute.correlation.samples <- function(params){
    L <- params$L
    iterations <- dim(L)[1]
    num.groups <- dim(L)[2]
    correlation <- array(dim=c(iterations, num.groups, num.groups))

    for(iter in 1:iterations){
      correlation[iter, , ] <- L[iter, , ] %*% t(L[iter, , ])
    }
    return(correlation)
}


#' Function to generate estimates of known group sizes
#' 
#' @export
#' @param params list of posterior samples of rho and L
#' @param population.names vector of population names
#' @param unknown.indices vector of indices of unknown groups
#' @param pop.col.name name of column in group.info containing population sizes
#' @param group.info data frame of group information
#' @return data frame of estimates of known group sizes
#' 
generate.spatial.estimates <- function(params, population.names, unknown.indices, pop.col.name, group.info){
  correlation <- compute.correlation.samples(params)
  m.corr <- apply(correlation, c(2, 3), mean)
  row.names(m.corr) <- population.names
  colnames(m.corr) <- population.names

  num.groups <- nrow(group.info)
  
  estimated.sizes <- matrix(0, nrow = num.groups, ncol = length(population.names))
  for(group in 1:nrow(estimated.sizes)){
    after.warmup.rho_j <- as.data.frame(params$rho_j[, group , ])
    population.size <- group.info[group, {{pop.col.name}}]
    known.population.names <- population.names[-{{unknown.indices}}]

    raw.known.sizes <- group.info[group, {{ known.population.names }} ]

    correlatedScaling(after.warmup.rho_j, correlation,
                      raw.known.sizes,
                      population.size) -> after.warmup.rho.scaled

    mean.rho <- colMeans(after.warmup.rho.scaled, na.rm=T)
    estimated <- round(exp(mean.rho)*population.size)

    estimated.sizes[group, ] <- estimated
  }

  estimated.sizes <- as.data.frame(estimated.sizes)
  colnames(estimated.sizes) <- population.names
  return(estimated.sizes)
}