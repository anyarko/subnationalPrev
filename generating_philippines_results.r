# Loading packages needed for modelling
# library(here)
library(rio)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)

options(scipen=999)

suppressPackageStartupMessages(source("utils.R"))

# to have access to the data folder and saved posterior samples
setwd("../")

# Loading the raw survey data
survey.org <- import("data/PH2205043401CP_CAM_Full_v1_2023Feb03.csv")


# selecting relevant data from the survey
survey <- survey.org %>% 
select(
region = RP.Reg,    
municipality = RP.City,
province = RP.Prov,   
age = S4a, 
work_mode = S5,
occupation = S5a,
montly_household_income = S7,
highest_education = S8,
household_size = S9num,
gender = S2,
doctors_known = `Grid_Q6a[{_1}].Q6a`,
lawyers_known = `Grid_Q6a[{_2}].Q6a`,
teachers_known = `Grid_Q6a[{_3}].Q6a`,
drivers_known = `Grid_Q6a[{_4}].Q6a`,
nurses_known = `Grid_Q6a[{_5}].Q6a`,
guards_known = `Grid_Q6a[{_6}].Q6a`,
students_known = `Grid_Q6a[{_7}].Q6a`,
constructions_known = `Grid_Q6a[{_8}].Q6a`,
bpos_known = `Grid_Q6a[{_9}].Q6a`,
sarisari_known = `Grid_Q6a[{_10}].Q6a`,
socials_known = `Grid_Q6a[{_11}].Q6a`,
aged_0_2_known = `Grid_Q6b[{_1}].Q6`,
aged_3_5_known = `Grid_Q6b[{_2}].Q6`,
aged_6_12_known = `Grid_Q6b[{_3}].Q6`,
aged_13_17_known = `Grid_Q6b[{_4}].Q6`,
perps_known = Q9_1,
live_perps_known = Q10_1,
children_known = `Grid_Q12[{_4}].Q12`,
live_children_known = `Grid_Q12[{_3}].Q12`)


# the data is sparse and people have a naturally cognitive limit on the number of 
# people they can know so limit the number of to below the 97.5 percentile for the
# occupation counts
survey.truncated <- survey %>% mutate_at(vars(doctors_known:aged_13_17_known), 
            funs(ifelse(. > quantile(., probs = seq(0,1, length.out=41))[39], 
                            quantile(., probs = seq(0,1, length.out=41))[39], .))) 


# calculating the average number of adults in households in each municipality from the survey data. 
# the final sizes are rounded to the nearest integer because integers are needed for modelling ARD
municipality.household.sizes <- survey %>% select(municipality, province, household_size) %>% 
    mutate(household_size = pmax(household_size, 1)) %>% 
    group_by(municipality, province) %>% summarise(mean.household.size = mean(household_size, na.rm = TRUE) %>% round) %>% as.data.frame


# dealing with people who reported victims and no perpetrators or vice versa and 
# the first combination of the non-live and live by just summing the two
survey.truncated.reporting <- survey.truncated
survey.truncated.reporting <- survey.truncated.reporting %>% 
                              mutate(combined_perps_1 = perps_known + live_perps_known, 
                                     combined_victims_1 = children_known + live_children_known)


# identifying respondents who reported any victims or perpetrators
any.victims <- which(survey.truncated.reporting$live_children_known != 0 | survey.truncated.reporting$children_known != 0)
any.perps <- which(survey.truncated.reporting$live_perps_known != 0 | survey.truncated.reporting$perps_known != 0)


# for respondents reporting perps and no victims impute the average number of victims per trafficker
for(respondent in setdiff(any.perps, any.victims)){
    respondent.municipality <- survey.truncated.reporting[respondent, "municipality"]
    respondent.province <- survey.truncated.reporting[respondent, "province"]

    survey.truncated.reporting[respondent, c("live_children_known", "children_known")] <- 4
    survey.truncated.reporting[respondent, c("combined_victims_1")] <- 4
}


# for respondents reporting victims and no perps impute the municipal average number of adults in households
for(respondent in setdiff(any.victims, any.perps)){
    respondent.municipality <- survey.truncated.reporting[respondent, "municipality"]
    respondent.province <- survey.truncated.reporting[respondent, "province"]

    municipality.mean.household.size <- municipality.household.sizes %>% 
        filter(municipality == respondent.municipality & province == respondent.province) %>% 
        pull(mean.household.size)

    survey.truncated.reporting[respondent, c("live_perps_known", "perps_known")] <- municipality.mean.household.size
    survey.truncated.reporting[respondent, c("combined_perps_1")] <- municipality.mean.household.size
}


# saving cleaned data
# survey.truncated.reporting %>% write.xlsx("ARD Truncated Reporting.xlsx")


# loading cleaned survey data
survey <- import("data/ARD Truncated Reporting.xlsx")


# loading municipality populations
population.info <- import("data/current_municipal_data.csv") %>% 
                   select(municipality=ADM3_EN, province=ADM2_EN, `Household Population`)


# grouping survey respondents by municipality. Province is required to distinguish municipalities
survey.pooling <- survey %>% group_by(municipality, province) %>%
                             mutate( group = cur_group_id() ) %>% 
                             mutate_at(vars(group), as.factor) %>% 
                             as.data.frame()


# appending municipality populations to survey data
pooling.group.names <- survey.pooling %>% select(municipality, province, group) %>% 
                       distinct() %>% left_join(population.info) %>% 
                       select(-municipality, -province) 

# municipalities that couldn"t be found automatically
pooling.group.names[1, "Household Population"] <- 801439
pooling.group.names[2, "Household Population"] <- 1659025
pooling.group.names[3, "Household Population"] <- 1837785
pooling.group.names[4, "Household Population"] <- 246743
pooling.group.names[55, "Household Population"] <- 115127
pooling.group.names[140, "Household Population"] <- 695410


survey.pooling <- survey.pooling %>% left_join(pooling.group.names, by=c("group"))


# This is the group number for each respondent representing the municipality they belong to
grouping <- survey.pooling %>% pull(group) %>% as.vector() %>% as.integer()


# maximum of live and non-live known perpetrators and victims is used since the individual
# categories are too sparse to model
df <- survey %>% 
      mutate(combined_perps_2 = pmax(perps_known, live_perps_known),
             combined_victims_2 = pmax(children_known, live_children_known), .keep="unused") %>% 
      select(drivers_known, constructions_known, sarisari_known, teachers_known,
            combined_perps_2, combined_victims_2) %>% 
      mutate_all(as.integer) %>% 
      as.matrix()


# population of the Philippines at time of survey
num.filipinos <- 114597229


# creating object for the names of the populations
population.names <- df %>% colnames() %>% gsub(pattern="_known", replacement = "")
known.population.names <- population.names[!grepl("combined", population.names)]


# loading the shapefile of the Philippines and indexing the municipalities surveyed from
# all the municipalities
library(sf)
adj <- st_read("data/AOI/current_municipalities.shp")
all.municipalities <- adj %>% st_drop_geometry() %>% select(ADM3_EN, ADM2_EN) %>% mutate(idx = row_number())
surveyed.municipalities <- survey.pooling %>% select(municipality, province) %>% distinct() %>% 
                           left_join(all.municipalities, c("municipality"="ADM3_EN", "province"="ADM2_EN"))

surveyed.municipalities[1, "idx"] <- 413
surveyed.municipalities[2, "idx"] <- 309
surveyed.municipalities[3, "idx"] <- 1
surveyed.municipalities[4, "idx"] <- 410
surveyed.municipalities[55, "idx"] <- 812

group.info <- surveyed.municipalities %>% cbind(pooling.group.names) %>% 
              mutate_at(vars(group), as.integer) %>% arrange(group)
group.info %>% head


# running the model. Can install the package and skip ahead to line 202 if saved 
# model posterior samples are available

# to build the subnationalPrev package
#  1. navigate to the directory with the "subnationalPrev" folder
#  2. run the following in an R console:
#  devtools::build("subnationalPrev")
#  this will generate a tar.gz file
#  3. install from terminal with admin rights the command:
#  RCMD INSTALL subnationalPrev_1.0.tar.gz
#  this gives access to the runSampling, posteriorRhoScaling, functions and the easy 
#  generateSpatialEstimates (without conf intervals)

library(subnationalPrev)
fit <- runSampling(ard=df, grouping=grouping, num.iterations=5000, warmup=3000, 
                   num.chains=2)                   

# loading the already fitted model and posterior samples
library(subnationalPrev)
load(file = "fits/neg_binomial_73.Rdata")
params <- rstan::extract(fit)


# the known groups were estimated from the occurrence in the households surveyed in philippines_occupation_estimation.r
known.group.sizes <- import("data/Estimated Known Group Sizes.xlsx") 
known.group.sizes[84, "group"] <- 1


# appending the size of the known groups to each surveyed municipality
group.info <- group.info %>% arrange(group) %>% cbind(known.group.sizes %>% select(drivers:teachers))
group.info %>% head


# function to get the known group sizes for each group.id used in the stan modelling
getKnownGroupSizes <- function(group.id){
  sizes <- known.group.sizes %>% filter(group == {{group.id}})
  sizes <- sizes[, known.population.names] %>% as.vector
  return(sizes)
}


# loop to generate estimates of all group sizes and 95% credible interval
estimated.sizes <- matrix(0, nrow = 149, ncol = 19)
for(group in 1:nrow(estimated.sizes)){
  print(paste("Currently scaling group", group))

  after.warmup.rho_j <- params$rho_j[, group , ] %>% as.data.frame * -1
  population.size <- known.group.sizes$"Household Population"[group]

  raw.known.sizes <- getKnownGroupSizes(group)  

  posteriorRhoScaling(after.warmup.rho_j,
                    raw.known.sizes,
                    population.size) -> after.warmup.rho.scaled

  mean.rho <- after.warmup.rho.scaled %>% colMeans(na.rm=T)
  estimated <- (exp(mean.rho)*population.size) %>% round 

  estimated.sizes[group, 1] <- group
  for(pop in 1:6){
    q <- quantile(after.warmup.rho.scaled[, pop], probs=c(0.025, 0.975), na.rm=T)
    estimated.sizes[group, (3*pop) - 1] <- (exp(q["2.5%"])*population.size) %>% round
    estimated.sizes[group, 3*pop] <- estimated[pop]
    estimated.sizes[group, (3*pop) + 1] <- (exp(q["97.5%"])*population.size) %>% round
  }
}

estimated.sizes <- estimated.sizes %>% as.data.frame 
colnames(estimated.sizes) <- c("group", rep(paste0("estimated_", population.names), each=3))

# adding the labels to the interval boundary columns
for(n in 1:length(colnames(estimated.sizes))){
  if(n == 1){
    next
  }
  if(n %% 3 == 1){
    colnames(estimated.sizes)[n] <- paste0(colnames(estimated.sizes)[n], "_95%_upper")
  }
  if(n %% 3 == 2){
    colnames(estimated.sizes)[n] <- paste0(colnames(estimated.sizes)[n], "_95%_lower")
  }
}


colnames(estimated.sizes)
library(openxlsx)
estimated.sizes %>% write.xlsx("results/Final Estimated Group Sizes Neg Binomial 5k w CI.xlsx", rowNames=F)
