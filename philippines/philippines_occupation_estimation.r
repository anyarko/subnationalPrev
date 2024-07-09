library(missMethods)
library(dplyr)
library(tidyr)
library(rio)
library(stringr)
library(openxlsx)

suppressPackageStartupMessages(source("utils.R"))

# to have access to the data folder
setwd("../")

# loading the raw survey data
survey.org <- import("data/PH2205043401CP_CAM_Full_v1_2023Feb03.csv")

# select the survey responses for household occupation information. The respondent is also
# counted here so no need to include their earlier response
survey <- survey.org %>% 
select(
municipality = RP.City,
province = RP.Prov,
ends_with("S9_WS"),
ends_with("S9_WSa")) %>%
mutate_at(vars(municipality), correctCityName) %>% 
mutate_at(vars(province), correctProvinceName)

# get the number of people in households for each municipality
occupation <- import("data/OccupationData.csv") %>% select(8:65) %>% impute_median

# identify the students 18 years and above from the responses
occupation.students <- survey %>% 
select(
municipality, province,
ends_with("S9_WS")) %>% 
mutate_at(vars(ends_with("S9_WS")), ~ ifelse(grepl("[24]+", .), "{_11}", "")) %>% 
pivot_longer(cols = ends_with("S9_WS"), names_to = "question", values_to = "occupation") %>%
select(-question)

# identify the other occupations from the responses and combine with the students
occupation.responses <- survey %>% select(
municipality, province,
ends_with("S9_WSa")) %>% 
pivot_longer(cols = ends_with("S9_WSa"), names_to = "question", values_to = "occupation") %>% 
select(-question) %>% rbind(occupation.students) %>%
filter(occupation != "") %>%
group_by(municipality, province, occupation) %>%
summarise(number = n()) %>%
arrange(municipality, province, occupation) %>% 
pivot_wider(names_from = occupation, values_from = number) 

# get the total number of people in households asked for each municipality
total.asked <- survey.org %>%
select(
municipality = RP.City,
province = RP.Prov,
S9num, 
S10num) %>% 
mutate_at(vars(municipality), correctCityName) %>% 
mutate_at(vars(province), correctProvinceName) %>% 
mutate(total_asked = S9num + S10num, .keep="unused") %>% 
group_by(municipality, province) %>% 
summarise_all(sum)

# the code to identify each known group
occupation.code <- c("{_1}" = "doctors", "{_2}" = "lawyers", "{_3}" = "teachers", "{_4}" = "drivers",
                     "{_5}" = "nurses", "{_6}" = "guards", "{_7}" = "constructions", "{_8}" = "bpos",
                     "{_9}" = "sarisari", "{_10}" = "other", "{_11}" = "students") 

colnames(occupation.responses)[3:13] <- occupation.code[colnames(occupation.responses)[3:13]] %>% unname

occupation.responses %>% as.data.frame() %>% head

# reporting no occupations as 0 instead of NA
occupation.responses <- occupation.responses %>% mutate_at(vars(-municipality, -province), ~ ifelse(is.na(.), 0, .)) %>% 
                       as.data.frame %>% left_join(total.asked)

# getting the estimated number of people in each known group based on percentage of responses
occupation.responses <- occupation.responses %>% 
mutate(available_known_groups = rowSums(across(drivers:lawyers) >0 )) 
occupation.responses %>% write.xlsx("Household Occupation Frequencies.xlsx", rowNames=F)
occupation.responses %>% head

population.info <- import("data/current_municipal_data.csv") %>% 
                   select(municipality=ADM3_EN, province=ADM2_EN, `Household Population`)

survey <- survey %>% select(municipality, province) %>% distinct() %>% 
          arrange(municipality, province) %>% 
          left_join(population.info) 

survey %>% filter(is.na(`Household Population`)) %>% select(municipality, province)

survey[which(survey$"Household Population"  %>% is.na),]

survey[35, "Household Population"] <- 1659025
survey[42, "Household Population"] <- 695410
survey[51, "Household Population"] <- 1837785
survey[53, "Household Population"] <- 246743
survey[55, "Household Population"] <- 801439
survey[84, "Household Population"] <- 115127

occupation.responses %>% left_join(survey) %>% 
mutate_at(vars(other:lawyers), ~ ((. / total_asked)* `Household Population`) %>% round) %>% 
write.xlsx("Estimated Known Group Sizes.xlsx", rowNames=F)
