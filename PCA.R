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
#   PCA plots
# 
# ******************************************************


## Load package(s):
library(tidyverse)
library(FactoMineR) # For PCA plots
library(ordr) # manipulating data objects

library(devtools) # to install ggord
install_github("fawda123/ggord") # install ggord from GitHub
library(ggord) # manipulating data objects




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

## join this aux. dataset to main PCA data
PCA_data <- left_join(PCA_data, ints_standard, by = "early_interval")

## Remove those that have NA (stratigraphic range is too long)
PCA_data <- na.omit(PCA_data, epoch)

## Filter to the 4 specific intervals
## For example, the early Cretaceous:
PCA_data_earlyK <- PCA_data %>% filter(epoch == "Early Cretaceous") 



# PCA: analysis and plots --------------------------------------------------

## Perform a separate PCA and plot for each interval and 



eK_pca <- PCA_data_earlyK[, 2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() #%>% # if error, check package 'ordr' is installed correctly
  #mutate_rows(group = PCA_data_eEC$family)

confellip_eK <- eK_pca %>% #na.omit(diet_group) %>% 
  ordr::ggbiplot(aes(color = family)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = family), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC")) + ## add more colours if >5 families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC")) ## add more colours if >5 families!
confellip_eK



## repeat for each interval!



