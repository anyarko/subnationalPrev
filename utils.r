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
               
generate.binomial.ard <- function(num.respondents, num.subpopulations, total.pop.size, p.k){
  rho_j <- rbinom(n=num.subpopulations, size=total.pop.size, prob = p.k)

  degree_mu <- runif(1, min=3, max=5)
  degree_sd <- runif(1, min=0.25, max=1.5)
  degrees <- rlnorm(num.respondents, meanlog=degree_mu, 
                    sdlog=degree_sd) %>% round

  degrees <- pmin(degrees, total.pop.size)

  ard <- matrix(0, nrow=num.respondents, ncol=num.subpopulations)
  for(respondent in 1:num.respondents){
    ard[respondent,] <- rbinom(num.subpopulations, size=degrees[respondent], 
                              prob=rho_j/total.pop.size)
  }
  return(ard)
}

generate.corr.nb.ard <- function(num.respondents, num.subpopulations, total.pop.size, p.k, w){
  rho_j <- log(p.k)
  delta <- rweibull(num.respondents, shape=35, scale=10)

  tau_n <- rep(0.2, num.subpopulations)
  mu  <- log(1 / sqrt(1 + (tau_n)**2))
  tau <- sqrt(log(1 + (tau_n)**2))

  corr.matrix <- matrix(c(1, 0.82, 0.35, 0.1, 0.01, 0.03, 
                          0.82, 1, 0.25, 0.1, 0.02, 0.05,
                          0.35, 0.25, 1, 0.3, 0.03, 0.06,
                          0.1, 0.1, 0.3, 1, 0.24, 0.27,
                          0.01, 0.02, 0.03, 0.24, 1, 0.9,
                          0.03, 0.05, 0.06, 0.27, 0.9, 1),
                          nrow=num.subpopulations, 
                          ncol = num.subpopulations)

  L <- t(chol(corr.matrix))

  rho_j <- matrix(rho_j, nrow=num.respondents, 
                          ncol=num.subpopulations, 
                          byrow = T)

  norm.errors <- matrix(rnorm(num.respondents * num.subpopulations, mean=0, sd=1), 
      ncol=num.subpopulations)

  lambda <- exp(rho_j - matrix(rep(delta, times=num.subpopulations), ncol=num.subpopulations) + 
              matrix(rep(mu, each=num.respondents), ncol=num.subpopulations) +
              t(diag(tau) %*% L %*% t(norm.errors))) * total.pop.size
  

  ard <- matrix(0, nrow=num.respondents, ncol=num.subpopulations)      
  for(respondent in 1:num.respondents){
    ard[respondent,] <- rnbinom(num.subpopulations, mu=lambda[respondent,], size=w)
  }
       
  return(ard)
}

generate.uncorr.nb.ard <- function(num.respondents, num.subpopulations, total.pop.size, p.k, w){
  rho_j <- log(p.k)

  delta <- rweibull(num.respondents, shape=16, scale=4)

  tau_n <- rep(0.2, num.subpopulations)
  mu  <- log(1 / sqrt(1 + (tau_n)**2))
  tau <- sqrt(log(1 + (tau_n)**2))
  tau  <- matrix(tau, nrow = num.respondents, ncol=num.subpopulations, byrow=T)

  rho_j <- matrix(rho_j, nrow=num.respondents, 
                          ncol=num.subpopulations, 
                          byrow = T)
  norm.errors <- matrix(rnorm(num.respondents * num.subpopulations, mean=0, sd=1), 
      ncol=num.subpopulations)

  lambda <- exp(rho_j + matrix(rep(delta, times=num.subpopulations), ncol=num.subpopulations) + 
              matrix(rep(mu, each=num.respondents), ncol=num.subpopulations) +
              (tau * norm.errors))        

  ard <- matrix(0, nrow=num.respondents, ncol=num.subpopulations)      
  for(respondent in 1:num.respondents){
    ard[respondent,] <- rnbinom(num.subpopulations, mu=lambda[respondent,], size=w)
  }

  browser()
  return(ard)
}

posteriorRhoScaling <- function(rho, raw_known_sizes,
                              population_size, correlation=NULL){
  iterations <- nrow(rho)
  K <- ncol(rho)
  new_rho <- rho
  known_ind <- which(raw_known_sizes > 0)
    
  for(k in 1:K){
    for(iter in 1:iterations){
      l_known_ind <- setdiff(known_ind, k)
      num_known <- length(l_known_ind)

      if(!is.null(correlation)){
        w <- correlation[iter, k, l_known_ind]
        w[w < 0] <- 0 
        w <- w * ( num_known / sum(w) )
      }
      
      known_sizes <- raw_known_sizes[-k]
      known_sizes <- known_sizes[which(known_sizes > 0)]
      
      if(!is.null(correlation)){
        C_m  <- log( (1/(num_known)) * 
                     sum( (exp(rho[iter, l_known_ind]) * w ) / (unlist(known_sizes)/population_size)) )
      }else{
        C_m  <- log( (1/(num_known)) * 
                     sum( exp(rho[iter, l_known_ind]) / (unlist(known_sizes)/population_size)) )
      }
      new_rho[iter, k] <- rho[iter, k] - C_m
    }
  }
  return(new_rho)
}

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

display.min.max.probs.gamma <- function(mu, sigma){
  shape.rho <- (mu**2)/(sigma**2)
  rate.rho <- mu/(sigma**2)
  print(exp(-rgamma(1000, shape.rho, rate.rho)) %>% min)
  print(exp(-rgamma(1000, shape.rho, rate.rho)) %>% max)
}

display.counts.negbin <- function(mu, size){
  print(table(rnbinom(10000, mu=mu, size=size)))
}
