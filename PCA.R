# ******************************************************
#
#   Pterosaur climate niche evolution
#
#   MSc thesis 2024
#
#   Lead: Carolin Grether (carolin.grether@fau.de)
#   Supervisor: Emma Dunne (emma.dunne@fau.de)
# ______________________________________________________
#
#   Raincloud plots & PCA
# 
# ******************************************************


## Load package(s):
library(tidyverse)
library(FactoMineR) # For PCA plots

library(devtools)
install_github("fawda123/ggord") ## For PCA plots
library(ggord)



# 1. Organise data -----------------------------------------------------------


## Import species/climate data if not already loaded:
#species_climate <- read_csv("Data/climate/species_climate.csv")
glimpse(species_climate)

## Take only the columns we need for the PCA:
PCA_data <- subset(species_climate, select = c(occurrence_no,
                                               MAT, seasonal_temp, 
                                               MAP, seasonal_precip, 
                                               #group, ## ***** THIS WILL BE YOUR GROUPING VARIABLE **** 
                                               early_interval
))


## Load auxillery dataset to standardise the time intervals
ints_standard <- read_csv("Data/Input/ints_standard.csv")

## rename 'early interval'
PCA_data <- rename(PCA_data, interval_std = early_interval)

## join aux. dataset to main PCA data
PCA_data <- left_join(PCA_data, ints_standard, by = "interval_std")

## Filter to the 4 specific intervals
## For example, the early Cretaceous:
PCA_data_EK <- PCA_data %>% filter(std_interval == "Early Cretaceous") 



# PCA: analysis and plots --------------------------------------------------

## Perform a separate PCA and plot for each interval and 



eEC_pca <- PCA_data_eEC[, 2:5] %>%
  prcomp(scale = TRUE) %>%
  as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_eEC$diet_group)

confellip_eEC <- eEC_pca %>% na.omit(diet_group) %>% 
  ggbiplot(aes(color = group)) +
  theme_bw() +
  geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC"))
confellip_eEC







