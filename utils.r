library(stringr)
library(rio)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)

options(scipen=999)


killworth <- function(ard, known.group.indices, known.group.sizes, total.pop.size) {
  num.respondents <- nrow(ard)
  degrees <- vector(length=num.respondents)

  for(respondent in 1:num.respondents){
    degrees[respondent] <- total.pop.size * ( sum(ard[respondent, known.group.indices]) / sum(known.group.sizes) )
  }

  copy.ard <- ard
  knew.none <- which(degrees == 0)
  if(length(knew.none) == num.respondents){
    return(rep(0, ncol(ard) - length(known.group.indices)))
  }

  if(length(knew.none) > 0){ 
    degrees <- degrees[-knew.none]
    copy.ard <- ard[-knew.none,] 
  }

  if(is.matrix(copy.ard)){
    key.respondents <- copy.ard[,-known.group.indices]
  }else{
    key.respondents <- copy.ard[-known.group.indices]
  }

  if(is.vector(key.respondents)){
    return(pmin( key.respondents / degrees, 1))
  }

  unknown.group.estimates <- colSums(key.respondents) / sum(degrees)

  return(unknown.group.estimates)
}
               

generate.nb.ard <- function(num.respondents, num.subpopulations, total.pop.size, p.k, w){
  rho_j <- log(p.k)
  delta <- rnorm(num.respondents, mean=5.5, sd=1)

  rho_j <- matrix(rho_j, nrow=num.respondents, 
                          ncol=num.subpopulations, 
                          byrow = T)

  lambda <- exp(delta + rho_j)

  ard <- matrix(0, nrow=num.respondents, ncol=num.subpopulations)      
  for(respondent in 1:num.respondents){
    ard[respondent,] <- rnbinom(num.subpopulations, mu=lambda[respondent,], size=w)
  }
  return(ard)
}


display.error.contributions <- function(all.p.k, num.regions, correlation, 
                                        known.group.indices, population){
  for(region in 1:num.regions){
    true.p.k <- all.p.k[region,]
    rho <- rho_j[,region,] * -1
    nb.estimate <- posteriorRhoScaling(rho, 
                      population*all.p.k[region, known.group.indices], 
                      population)               
    mean.estimate <- colMeans(nb.estimate, na.rm=T) %>% pmin(0)                  

    region.error  <- mean(abs((exp(mean.estimate[-known.group.indices]) - true.p.k[-known.group.indices])) 
      / true.p.k[-known.group.indices])

    print(paste('Mean Absolute Relative Error for region', region, 'is', region.error))
  }
}


display.simulated.rho_j <- function(region, all.p.k, fit){
  parameter.names <- dimnames(fit)$parameters
  pars <- parameter.names[grepl(paste0('^rho_j\\[', region, ','), parameter.names) ]
  print(fit, pars=pars)
  print('True prevalence the region is: ')
  print(all.p.k[region,] %>% log())
}
