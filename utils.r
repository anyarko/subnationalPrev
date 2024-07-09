library(lubridate)
library(stringr)
library(purrr)
library(rio)
library(sf)
library(tibble)
library(dplyr)
library(tidyr)
library(reshape2)

options(scipen=999)

interchangeCityName <- function(df1, df2, col.name1, col.name2, full=FALSE){
    # Adjust the names of municipalities in df1 to match df2 if the
    # difference is just the change of ... City to City of ...
    municipalities.df1 <- df1 %>% pull(col.name1) %>% unique %>% sort
    municipalities.df2 <- df2 %>% pull(col.name2) %>% unique %>% sort
    missing.municipalities.df1 <- setdiff(municipalities.df1, municipalities.df2)
    missing.municipalities.df2 <- setdiff(municipalities.df2, municipalities.df1)
    if(full){ missing.municipalities.df2  <- municipalities.df2 }

    for(municipality in missing.municipalities.df1){
        if(length(grep(' City', municipality))){
            check <- paste('City of', sub(' City', '', municipality))
            if(check %in% missing.municipalities.df2){ 
              if(inherits(df1, 'sf')){
                df1[(st_drop_geometry(df1[, col.name1]) == municipality ), col.name1]  <- check
              }else{
                df1[(df1[, col.name1] == municipality ), col.name1]  <- check
              }
            }
        }
    }
    return(df1)
}


addCityOf <- function(df1, df2, col.name1, col.name2, full=FALSE){
    # Adjust the names of municipalities in df1 to match df2 if the
    # difference is just the change of ... City to City of ...
    municipalities.df1 <- df1 %>% pull(col.name1) %>% unique %>% sort
    municipalities.df2 <- df2 %>% pull(col.name2) %>% unique %>% sort
    missing.municipalities.df1 <- setdiff(municipalities.df1, municipalities.df2)
    missing.municipalities.df2 <- setdiff(municipalities.df2, municipalities.df1)
    if(full){ missing.municipalities.df2  <- municipalities.df2 }

    for(municipality in missing.municipalities.df1){
        if(length(grep('City', municipality)) == 0){
            check <- paste('City of', municipality)
            if(check %in% missing.municipalities.df2){ 
              if(inherits(df1, 'sf')){
                df1[(st_drop_geometry(df1[, col.name1]) == municipality ), col.name1]  <- check
              }else{
                df1[(df1[, col.name1] == municipality ), col.name1]  <- check
              }
            }
        }
    }
    return(df1)
}



wrong.names <- c('Tagoloan Ii'='Tagoloan II', "Brooke's Point"="Brooke's Point", 'Gen. S.k. Pendatun'='Gen. S.K. Pendatun',
                 'City of Ozamiz'='City of Ozamis', 'Pinamungajan'='Pinamungahan', 'Mataas Na Kahoy'='Mataasnakahoy', 
                 'City of Calbayog'='Calbayog City', 'Cordoba'='Cordova', 'Jetafe'='Getafe', 'Pozzorubio'='Pozorrubio',
                 'Manila'='City of Manila', 'Amai Manabilang'='Bumbaran', 'Alfonso Castañeda'='Alfonso Castaneda',
                 'Bagbag'='Bagabag', 'Barauen'='Burauen')

manila.districts <- c('Ermita', 'Intramuros', 'Binondo', 'Malate', 'Paco', 'Pandacan', 'Port Area', 'San Juan',
                      'Quiapo', 'Sermita', 'Tondo I/Ii', 'Tondo', 'Metro Manila', 'Santa Mesa', 'Sampaloc')

correctCityName <- function(municipality.name){
  if(is.na(municipality.name)){ return(NA)}
  municipality.name <- municipality.name %>% str_to_title
  if(length(grep('Of', municipality.name))){
    municipality.name <- sub('Of', 'of', municipality.name)
  }
  if(length(grep(' De', municipality.name))){
    municipality.name <- sub(' De', ' de', municipality.name)
  }
  if(length(grep('\\.\\.', municipality.name))){
    municipality.name <- gsub('\\.\\.', '', municipality.name)
  }
  if(length(grep('[\n\r]+', municipality.name))){
    municipality.name <- gsub('[\n\r]+', ' ', municipality.name)
  }
  if(length(grep('\\([A-Za-z &-ñ.)(]+', municipality.name))){
    municipality.name <- sub('\\([A-Za-z &-ñ.)(]+', '', municipality.name)
  }
  municipality.name <- municipality.name %>% str_squish
  if(municipality.name %in% names(wrong.names)){
    return(wrong.names[municipality.name])
  }
  if(municipality.name %in% manila.districts){
    return("City of Manila")
  }
  num <- regexpr('-', municipality.name)
  if(num > 0){
    substr(municipality.name, num+1, num+1) <- substr(municipality.name, num+1, num+1) %>% str_to_upper()
  }
  num <- regexpr('\'', municipality.name)
  if(num > 0){
    substr(municipality.name, num+1, num+1) <- substr(municipality.name, num+1, num+1) %>% str_to_upper()
  }
  return(municipality.name)
}
correctCityName <- Vectorize(correctCityName)


correctProvinceName <- function(province.name){
  if(province.name %>% is.na){ return(NA)}
  province.name <- province.name %>% str_to_title()

  if(length(grep(' Del', province.name))){
    province.name <- sub(' Del', ' del', province.name)
  }
  if( province.name == 'Mountain Province'){
    return(province.name)
  }
  if(length(grep('Province', province.name))){
    province.name <- sub('Province', '', province.name)
  }
  if(length(grep('\r\n', province.name))){
    province.name <- sub('\r\n', '', province.name)
  }
  if(length(grep(' \\([A-Za-zñ -.)(]+', province.name))){
    province.name <- sub(' \\([A-Za-zñ -.)(]+', '', province.name)
  }
  if(substr(province.name, 1, 3) == 'Ncr'){
    substr(province.name, 1, 3) <- 'NCR'
  }
  if(length(grep('Of', province.name))){
    province.name <- sub('Of', 'of', province.name)
  }
  if(province.name == "City of Cotabato"){
    return("Cotabato City")
  }
  if(province.name == "North Cotabato"){
    return("Cotabato")
  }
  if(province.name == "Western Samar"){
    return("Samar")
  }
  if(province.name == "Tawi-tawi"){
    return("Tawi-Tawi")
  }
  if(province.name %in% c("NCR, First District", 'Manila') || startsWith(province.name, "NCR, City of Manila")){
    return("NCR, City of Manila, First District")
  }
  if(province.name == "Davao De Oro"){
    return("Compostela Valley")
  }
  
  return(province.name)
}
correctProvinceName <- Vectorize(correctProvinceName)


region.names  <- c("NCR"="National Capital Region", "CAR"="Cordillera Administrative Region", "REGION I"="Region I", "REGION II"="Region II", 
                   "REGION III"="Region III", "REGION IV-A"="Region IV-A", "REGION IV-B"="Region IV-B", 
                   "REGION V"="Region V", "REGION VI"="Region VI", "REGION VII"="Region VII", 
                   "REGION VIII"="Region VIII", "REGION IX"="Region IX", "REGION X"="Region X",
                   "REGION XI"="Region XI", "REGION XII"="Region XII", "CARAGA"="Region XIII", 
                   "BARMM"="Autonomous Region in Muslim Mindanao")


correctRegionName <- function(region.name){
  if(region.name %in% names(region.names)){
    return(region.names[region.name])
  }
  return(region.name)
}
correctRegionName <- Vectorize(correctRegionName)


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
