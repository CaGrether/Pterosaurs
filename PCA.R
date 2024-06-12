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
#   PCA plots
# 
# ******************************************************


## Load package(s):
library(tidyverse)
library(FactoMineR) # For PCA plots

library(devtools)
install_github("fawda123/ggord") ## For PCA plots
library(ggord)
library(ordr) # manipulating data objects

#library(devtools) # to install ggord
#install_github("fawda123/ggord") # install ggord from GitHub
#library(ggord) # manipulating data objects




# 1. Organise data -----------------------------------------------------------


## Import species/climate data if not already loaded:
species_climate <- read_csv("Data/climate/species_climate.csv")
glimpse(species_climate)

ptero_grouping <- read.csv2("Data/Input/ptero_groups_copy.csv")

## merge species climate and group data here
# only take species that we're using
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

unique(species_climate_group$family)

## Take only the columns we need for the PCA:
PCA_data <- subset(species_climate_group, select = c(occurrence_no,
                                               MAT, seasonal_temp, 
                                               MAP, seasonal_precip, 
                                               family, ## ***** THIS WILL BE YOUR GROUPING VARIABLE **** 
                                               early_interval
))


## Load auxillery dataset to standardise the time intervals
ints_standard <- read_csv2("Data/Input/ints_standard_copy.csv") # manually fixed data

## rename 'early interval'
PCA_data <- rename(PCA_data, interval_std = early_interval)
ints_standard <- rename(ints_standard, interval_std = early_interval) # added this

## join aux. dataset to main PCA data
PCA_data <- left_join(PCA_data, ints_standard, by = "interval_std") # only works with line above

## Filter to the 4 specific intervals
## For example, the early Cretaceous:
#PCA_data_EK <- PCA_data %>% filter(interval_std == "Early Cretaceous") #changed std_interval to interval_std
## join this aux. dataset to main PCA data
#PCA_data <- left_join(PCA_data, ints_standard, by = "interval_std") ??????

## Remove those that have NA (stratigraphic range is too long)
PCA_data <- na.omit(PCA_data, epoch)

## add empty row for numeric family values
PCA_data$family_num<- rep(NA,length(PCA_data$family))

## Filter to the 4 specific intervals
## For example, the early Cretaceous:
PCA_data_earlyK <- PCA_data %>% filter(epoch == "Early Cretaceous") 
PCA_data_earlyJ <- PCA_data %>% filter(epoch == "Early Jurassic") 
PCA_data_midJ <- PCA_data %>% filter(epoch == "Middle Jurassic")
PCA_data_lateT <- PCA_data %>% filter(epoch == "Late Triassic") 

unique(PCA_data_earlyK$family) # check all groups, best is n<6

## make family groups numeric
# late T
PCA_data_lateT$family_num[which(PCA_data_lateT$family=="Pterosauria")] <- 1
PCA_data_lateT$family_num[which(PCA_data_lateT$family=="Campylognathoididae")] <- 2
PCA_data_lateT$family_num[which(PCA_data_lateT$family=="Lagerpetidae")] <- 3
# early J
PCA_data_earlyJ$family_num[which(PCA_data_earlyJ$family=="Macronychoptera")] <- 4
PCA_data_earlyJ$family_num[which(PCA_data_earlyJ$family=="Campylognathoididae")] <- 2
PCA_data_earlyJ$family_num[which(PCA_data_earlyJ$family=="Dimorphodontidae")] <- 5
PCA_data_earlyJ$family_num[which(PCA_data_earlyJ$family=="Rhamphorhynchidae")] <- 6
PCA_data_earlyJ$family_num[which(PCA_data_earlyJ$family=="Rhamphorhynchoidea")] <- 7
# middle J
PCA_data_midJ$family_num[which(PCA_data_midJ$family=="Rhamphorhynchidae")] <- 6
PCA_data_midJ$family_num[which(PCA_data_midJ$family=="Macronychoptera")] <- 4
PCA_data_midJ$family_num[which(PCA_data_midJ$family=="Anurognathidae")] <- 8
PCA_data_midJ$family_num[which(PCA_data_midJ$family=="Scaphognathidae")] <- 9
PCA_data_midJ$family_num[which(PCA_data_midJ$family=="Wukongopteridae")] <- 10
# early K
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Ornithocheiroidea")] <- 11
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Lophocratia")] <- 12
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Azhdarchoidea")] <- 13
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Boreopteridae")] <- 14
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Anurognathidae")] <- 15
PCA_data_earlyK$family_num[which(PCA_data_earlyK$family=="Rhamphorhynchidae")] <- 6

length(PCA_data_lateT)
# PCA: analysis and plots --------------------------------------------------

## Perform a separate PCA and plot for each interval and 



#eEC_pca <- PCA_data_eEC[, 2:5] %>%
  #prcomp(scale = TRUE) %>%
  #as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  #mutate_rows(group = PCA_data_eEC$diet_group)

#confellip_eEC <- eEC_pca %>% na.omit(diet_group) %>% 
  #ggbiplot(aes(color = group)) +
  #theme_bw() +
  #geom_rows_point() +
  #geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  #geom_cols_vector(color = "#444444") + # adds the arrows
  #scale_colour_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC")) +
  #scale_fill_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", "#572AAC"))
#confellip_eEC


lT_pca <- PCA_data_lateT[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() #%>% # if error, check package 'ordr' is installed correctly
#mutate_rows(group = PCA_data_eEC$family)

family_num <- PCA_data_lateT$family_num

confellip_lT <- lT_pca %>% #na.omit(diet_group) %>% 
  ordr::ggbiplot(aes(color = family_num)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = family_num), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#1E44AA", "#248528", "#FFA93D")) + ## add more colours if >5 families!
  scale_fill_manual(values = c( "#1E44AA", "#248528", "#FFA93D")) ## add more colours if >5 families!
confellip_lT



eK_pca <- PCA_data_earlyK[, c(2:5,10)] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() #%>% # if error, check package 'ordr' is installed correctly
#mutate_rows(group = PCA_data_eEC$family)

confellip_eK <- eK_pca %>% #na.omit(diet_group) %>% 
  ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = family_num)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = family_num), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", 
                                 "#572AAC","#D7E05A")) + ## add more colours if >5 families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA", "#248528", "#FFA93D", 
                               "#572AAC","#D7E05A")) ## add more colours if >5 families!
confellip_eK



## repeat for each interval!



