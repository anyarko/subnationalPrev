# loading_libs_compiling_model
library(dplyr)
library(rstan)
library(rio)
options(scipen=999)

source('utils.r')
model <- stan_model('neg_binomial_partial_pooling.stan')

known.group.indices <- c(1, 2, 3, 4)
w <- c(20, 20, 20, 20, 0.1, 0.1)
mu.rho <- c(2.5, 3.5, 4.5, 5, 5.5, 6.5)
num.subpopulations <- length(mu.rho)
sigma.rho <- rep(1, num.subpopulations)

shape.rho <- (mu.rho**2) / (sigma.rho **2)
rate.rho <- (mu.rho) / (sigma.rho **2)

unknown.group.indices <- setdiff(1:num.subpopulations, known.group.indices)
sim.num.respondents <- c(100, 250, 500)
population <- 1000000
num.regions <- 150

# all.seeds <- c(971, 828, 621)

# generate_group_size
# set.seed(871) on Albert's work assigned laptop
# rho_j <- matrix(0, ncol=num.subpopulations, nrow=num.regions)
# for(region in 1:num.regions){
#   rho_j[region,]  <- -rgamma(num.subpopulations, shape=shape.rho, rate=rate.rho)                                 
# }
# save(rho_j, file = 'rho_j.Rdata')
load('rho_j.Rdata')

all.p.k <- rho_j %>% exp() 
colnames(all.p.k) <- paste0('subpopulation_', 1:num.subpopulations)
all.p.k %>% head(20)

for(s in 1:20){
all.seeds <- c(871, 728, 521) - s
# generate_multiple_sets_of_regions_simulate
set.seed(all.seeds[2])
num.regions <- 150

spatial.results <- matrix(0, ncol=3) %>% as.data.frame()
colnames(spatial.results) <- c('num_respondents', 'mean_error', 'model')

for(num.respondents in sim.num.respondents){
  cat(paste('Currently simulating with', num.respondents, 'respondents'), '\n')
  # Generate NB data for multiple regions and estimate the prevalence
  ard <- matrix(0, ncol=num.subpopulations)

  for(region in 1:num.regions){
    p.k <- as.numeric(all.p.k[region, ])
    nb.ard <- generate.uncorr.nb.ard(num.respondents = num.respondents,
                              num.subpopulations = num.subpopulations,
                              total.pop.size = population,
                              p.k = p.k, 
                              w = w)   
    ard <- rbind(ard, nb.ard)
  }
  ard <- ard[-1,]

  data <- list(N = nrow(ard),
               K = ncol(ard),
               y = ard,
               J = num.regions, 
               jj = rep(1:num.regions, each=num.respondents))

  fit <- sampling(model, data = data, warmup=200,
          iter = 300, cores=20,
          chains = 2, seed=all.seeds[3], verbose=F, show_messages=F)

  params <- extract(fit)
  rho_j <- params$rho_j

  for(region in 1:num.regions){
    true.p.k <- all.p.k[region,]
    rho <- rho_j[,region,] * -1

    nb.estimate <- posteriorRhoScaling(rho, 
                      population*all.p.k[region, known.group.indices], 
                      population)               
    mean.estimate <- colMeans(nb.estimate, na.rm=T) %>% pmin(0)                  

    region.error  <- mean(abs(exp(mean.estimate[-known.group.indices]) - true.p.k[-known.group.indices]) 
      / true.p.k[-known.group.indices])
    
    spatial.results <- rbind(spatial.results, 
                            c(num.respondents, 
                              region.error,
                              'Negative_Binomial'))
  }

  for(region in 1:num.regions){
    true.p.k <- all.p.k[region,]
    res <- c(1:num.respondents + (region-1)*num.respondents)
    binomial.estimate <- killworth(ard[res, ], 
                known.group.indices=known.group.indices, 
                known.group.sizes=p.k[known.group.indices]*population,
                total.pop.size=population)

    region.error  <- mean(abs(binomial.estimate - true.p.k[-known.group.indices])/ true.p.k[-known.group.indices])
    spatial.results <- rbind(spatial.results, 
                            c(num.respondents, 
                              region.error,
                              'Killworth_MLE'))
  }                                
}

spatial.results <- spatial.results[-1,] %>% mutate_at(vars(-model), as.numeric)

write.csv(spatial.results, paste('seeds_',
    paste(all.seeds, collapse='_'),
    '_spatial_results.csv', sep=''), row.names=F)
}
