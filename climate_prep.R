# ******************************************************
#
#   Pterosaur climate niche evolution
#
#   MSc thesis 2025
#
#   Lead: Carolin M. Grether (carolin.grether@fau.de)
#   Supervisor: Emma Dunne (emma.dunne@fau.de)
# ______________________________________________________
#
#   Prepping the climate data for analyses
# 
# ******************************************************


## Load packages
library(tidyverse)  # data manipulation and organisation
library(desk) # for Jarque Bera test


## Climate data - Scotese model output
climate_data_raw <- read.csv("Data/climate/tfks_ma.csv")

## First let's remove the collections that climate data could not be obtained for:
climate_data_complete <- climate_data_raw %>% drop_na(jan_temp_mm_srf)
glimpse(climate_data_complete) # check


#### Calculate climate variables
climate_data_varis <- climate_data_complete %>%
  rowwise() %>%
  ## Mean annual temperature (MAT)
  mutate(MAT = mean(c_across(c(jan_temp_mm_srf:dec_temp_mm_srf)), na.rm = TRUE)) %>%
  ## Seasonal variation in temperature
  mutate(temp_max = max(c_across(c(jan_temp_mm_srf:dec_temp_mm_srf)))) %>% # month with max temp
  mutate(temp_min = min(c_across(c(jan_temp_mm_srf:dec_temp_mm_srf)))) %>% # month with min temp
  mutate(seasonal_temp = temp_max - temp_min) %>% # difference between max and min temp month
  ## Mean annual precipitation (MAP)
  mutate(MAP = mean(c_across(c(jan_precip_mm_srf:dec_precip_mm_srf)), na.rm = TRUE)) %>%
  ## Seasonal variation in precipitation
  mutate(precip_max = max(c_across(c(jan_precip_mm_srf:dec_precip_mm_srf)))) %>% # month with max precip
  mutate(precip_min = min(c_across(c(jan_precip_mm_srf:dec_precip_mm_srf)))) %>% # month with min precip
  mutate(seasonal_precip = precip_max - precip_min) # difference between max and min precip month

glimpse(climate_data_varis) # check they have been added as new columns


## Clean up the tibble so that it only contains the variables we need
climate_data <- climate_data_varis %>% select(collection_no, plng, plat, 
                                              MAT, seasonal_temp, 
                                              MAP, seasonal_precip)


## Import the cleaned data (species and body fossils only)
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20) # fix file path if necessary
glimpse(occurrences)

## Combine both of these data by matching the species names (accepted_name)
species_info <- select(occurrences, collection_no, occurrence_no, accepted_name, 
                       early_interval, late_interval, max_ma, min_ma)
glimpse(species_info)

## Now combine this with the climate data from above
species_climate <- left_join(species_info, climate_data, by = "collection_no")
View(species_climate)

# ## Remove the rows that have any NAs - these are collections that we don't have species data for
# species_climate <- species_climate_full %>% drop_na(accepted_name)


## Save a copy as a .csv:
write_csv(species_climate, "Data/climate/species_climate.csv")

## test for normality
# which rows to remove in seasonal_temp (contain NA)
which(is.na(species_climate$seasonal_temp))

# Jarque Bera test for normal distribution
jbMAT <- jb.test(species_climate$MAT)
jbMAP <- jb.test(species_climate$MAP)
jbT <- jb.test(species_climate[-c(1090,1123),]$seasonal_temp) 
jbP <- jb.test(species_climate$seasonal_precip)


