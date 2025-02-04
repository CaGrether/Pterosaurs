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
#   Climate ranges - Raincloud plots
# 
# ******************************************************

## Load package(s):
library(tidyverse)
library(ggplot2)
library(ggdist)


# 1. Organise data -----------------------------------------------------------

## Import species/climate data if not already loaded:
  species_climate <- read.csv("Data/climate/species_climate.csv")
  glimpse(species_climate)

## Pterosaur clades
  ptero_grouping <- read.csv2("Data/Input/azhd_and_pteran.csv") # Azhdarchoidea and Pteranodontoidea
  glimpse(ptero_grouping)

## merge species climate and group data here
  species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
  species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

  #unique(species_climate_group$family)


  cloud_data <- subset(species_climate_group, select = c(occurrence_no,   
                                                     MAT, seasonal_temp,
                                                     MAP, seasonal_precip,
                                                     two.groups, 
                                                     early_interval,
                                                     accepted_name
  ))

## Load auxillary dataset to standardise the time intervals
  ints_standard <- read_csv2("Data/Input/ints_standard_copy.csv") # manually fixed data

## rename 'early interval'
  cloud_data <- rename(cloud_data, interval_std = early_interval)
  ints_standard <- rename(ints_standard, interval_std = early_interval)

## join aux. dataset to main PCA data
  cloud_data <- left_join(cloud_data, ints_standard, by = "interval_std")

## Remove those that have NA (stratigraphic range is too long)
  cloud_data <- na.omit(cloud_data, epoch)
  glimpse(cloud_data)
  
  # rename data
  colnames(cloud_data)[6] <- "Pterosaur_taxa"
  colnames(cloud_data)
  
  median(cloud_data$MAT[which(cloud_data$Pterosaur_taxa=="Azhdarchoidea")])
  median(cloud_data$MAT[which(cloud_data$Pterosaur_taxa=="Pteranodontoidea")])
  
## Filter to the specific intervals (EK and LK)
  cloud_EK <- cloud_data %>% filter(epoch == "Early Cretaceous") 
  cloud_LK <- cloud_data %>% filter(epoch == "Late Cretaceous") 


# Raincloud plots -------------------------------------------------------------
## MAT EK
ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = MAT, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Mean Annual Temperature (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

## MAT LK
ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = MAT, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Mean Annual Temperature (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


# Mann-Whitney U test
Wilcox_MAT <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_MAT_EK <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_MAT_LK <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## seasonal temp EK
ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = seasonal_temp, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Seasonal Temperature (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

# seasonal temp LK
ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = seasonal_temp, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Seasonal Temperature (°C)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

# Mann-Whitney U test
Wilcox_T <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_T_EK <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_T_LK <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## MAP EK
ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = MAP, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Mean Annual Precipitation (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

## MAP LK
ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = MAP, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Mean Annual Precipitation (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


# Mann-Whitney U test
Wilcox_MAP <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_MAP_EK <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_MAP_LK <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## seasonal Precip EK
ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = seasonal_precip, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Seasonal Precipitation (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

## seasonal Precip LK
ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = seasonal_precip, fill = Pterosaur_taxa)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA,
    alpha = 0.6) + 
  geom_boxplot(
    aes(fill = Pterosaur_taxa, colour = Pterosaur_taxa),
    width = .25, 
    outlier.shape = NA,
    alpha = 0.15
  ) +
  geom_point(
    aes(colour = Pterosaur_taxa),
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  ) + 
  scale_color_manual(values = c("#FFA93D", "#572AAC")) +
  scale_fill_manual(values = c("#FFA93D", "#572AAC")) +
  labs(x = NULL, y = "Seasonal Precipitation (mm/day)") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


# Mann-Whitney U test
Wilcox_P <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_P_EK <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_P_LK <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)


###########################

### fix seasonal temp that is above 40°C
## why is their seasonal temp so high? Is there a bias?

# dataset to get info of raincloud taxa
fix_data <- occurrences_sp # rename to keep original
fix_data$accepted_name <- gsub(" ", "_", fix_data$accepted_name)
fix_data$mid_ma <- (fix_data$max_ma + fix_data$min_ma)/2 # add mean age

# merge with all raincloud taxa
occ_info <- fix_data %>% 
  select(accepted_name, occurrence_no, mid_ma, collection_name, paleolat, paleolng, 
         taxon_environment, lagerstatten, collection_no) # columns that are of interest
cloud_all <- cloud_data %>% 
  merge(., occ_info, by = "occurrence_no") %>% 
  select(accepted_name, everything()) %>% 
  arrange(accepted_name) 

# subset for only >40°C taxa
cloud_40 <- cloud_all %>% 
  subset(.,.$seasonal_temp>40)

####1. "Problematic" localities
# check unique localities - how many occ nos are in there and are they the same as the >40 taxa?
# length(which(cloud_40$collection_name == "Chaoyang pterosaurs (PROXY)"))
# length(which(cloud_all$collection_name == "Chaoyang pterosaurs (PROXY)")) # same length

## Chaoyang pterosaurs (PROXY)
cloud_all$occurrence_no[which(cloud_all$collection_name == "Chaoyang pterosaurs (PROXY)")] %in% cloud_40$occurrence_no
# all Chaoyang pterosaurs (PROXY) are >40

## Yuanjiawa, Dapingfang
cloud_all$occurrence_no[which(cloud_all$collection_name == "Yuanjiawa, Dapingfang")] %in% cloud_40$occurrence_no
# all are > 40

#### 2. Lagerstatten
cloud_all$occurrence_no[which(cloud_all$lagerstatten == "conservation")] %in% cloud_40$occurrence_no
# not all conservation Lagerstatten are >40

#### 3. cold MAT's 
cloud_all$occurrence_no[which(cloud_all$MAT < 10)] %in% cloud_40$occurrence_no
# not all MAT < 10 are in seasonal temp > 40