#' Function to generate estimates of known municipality sizes
#' 
#' @export
#' @param params list of posterior samples of rho and L
#' @param population.names vector of population names same as column names of municipality.info
#' @param unknown.indices vector of indices of unknown municipalitys
#' @param pop.col.name name of column in municipality.info containing population sizes
#' @param municipality.info data frame of municipality information in the sample order as the ard passed to run.sampling
#' @return data frame of estimates of group sizes for each municipality
#' 
generateSpatialEstimates <- function(params, population.names, unknown.indices, pop.col.name, municipality.info){
  num.municipalities <- nrow(municipality.info)
  
  estimated.sizes <- matrix(0, nrow = num.municipalities, ncol = length(population.names))
  for(municipality in 1:num.municipalities){
    after.warmup.rho_j <- as.data.frame(params$rho_j[, municipality , ]) * -1
    population.size <- municipality.info[municipality, {{pop.col.name}}]
    known.population.names <- population.names[-{{unknown.indices}}]

    raw.known.sizes <- municipality.info[municipality, {{ known.population.names }} ]

    posteriorRhoScaling(after.warmup.rho_j,
                      raw.known.sizes,
                      population.size) -> after.warmup.rho.scaled

    mean.rho <- colMeans(after.warmup.rho.scaled, na.rm=T)
    estimated <- round(exp(mean.rho)*population.size)

    estimated.sizes[municipality, ] <- estimated
  }

  estimated.sizes <- as.data.frame(estimated.sizes)
  colnames(estimated.sizes) <- population.names
  return(estimated.sizes)
}