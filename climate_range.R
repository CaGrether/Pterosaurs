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
library(ggpubr) # for plotting

library(geoscale) # for plotting with the geological time scale on the x-axis (uses base R syntax)
library(rgplates) # paleogeographic reconstructions

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

## get collection names
  # occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20)
  # occurrences_sp <- occurrences %>% filter(accepted_rank == "species")
  # 
  # forfourty_sp <- occurrences_sp %>% 
  #   select(collection_name, occurrence_no, lat, lng, paleolat, paleolng, early_interval, late_interval, max_ma, min_ma) %>% 
  #   distinct(collection_name, .keep_all = TRUE) %>% 
  #   na.omit(collection_name)

## data for raincloud plots
  cloud_data <- subset(species_climate_group, select = c(occurrence_no,   
                                                     MAT, seasonal_temp,
                                                     MAP, seasonal_precip,
                                                     two.groups,
                                                     collection_no,
                                                     early_interval,
                                                     accepted_name,
                                                     plat, plng
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
  
  # find out medians
  median(cloud_data$MAT[which(cloud_data$Pterosaur_taxa=="Azhdarchoidea")])
  median(cloud_data$MAT[which(cloud_data$Pterosaur_taxa=="Pteranodontoidea")])
  
## Filter to the specific intervals (EK and LK)
  cloud_EK <- cloud_data %>% filter(epoch == "Early Cretaceous") 
  cloud_LK <- cloud_data %>% filter(epoch == "Late Cretaceous") 


# Raincloud plots -------------------------------------------------------------
## MAT EK
mat_EK <- ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = MAT, fill = Pterosaur_taxa)) + 
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
mat_EK

## MAT LK
mat_LK <- ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = MAT, fill = Pterosaur_taxa)) + 
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
mat_LK

# Mann-Whitney U test
Wilcox_MAT <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_MAT_EK <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_MAT_LK <- wilcox.test(MAT ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## seasonal temp EK
t_EK <- ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = seasonal_temp, fill = Pterosaur_taxa)) + 
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
t_EK

## find location of problematic species T > 40°C

# ----------------- method without corrections (species in ocean)
# subset species
# hot_sp <- cloud_EK[which(cloud_EK$seasonal_temp>40),]
# # assign EK age
# hot_sp$age <- 121
# # reconstruct model
# palgeoEK_hot <- reconstruct("plate_polygons", age = 121, model="MERDITH2021")
# ## map theme
# palaeomap_theme <- theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
#                                            axis.title.y=element_blank(), axis.text.y=element_blank(),
#                                            axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
#                                            legend.title=element_blank())
# Map_hot <-  ggplot() +
#   ## Landmasses
#   geom_sf(data = palgeoEK_hot, colour = "grey75", fill = "grey75") +
#   ## occurrence data
#   geom_point(data = cloud_EK, aes(x = plng, y = plat), color = "#FF8899", size = 4,  alpha = 0.8) + 
#   ## title 
#   ggtitle("Hot species") +
#   ## theme
#   palaeomap_theme
# Map_hot

# ------------------ corrections applied

## code to fix locations
# Data
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20)
occurrences_sp <- occurrences %>% filter(accepted_rank == "species")

## take necessary columns
taxon_inf <- select(occurrences_sp, collection_no, occurrence_no, accepted_name, 
                    early_interval, late_interval, max_ma, min_ma, formation, lng, lat)
# rename object
colls <- taxon_inf

# add mean age to each collection
colls$mid_ma <- (colls$max_ma + colls$min_ma)/2

# add strat info to filter for EK
## rename 'early interval'
colls <- rename(colls, interval_std = early_interval)

## join aux. dataset to main PCA data
colls <- left_join(colls, ints_standard, by = "interval_std")

## Remove those that have NA (stratigraphic range is too long)
colls <- na.omit(colls, epoch)
glimpse(colls)

## Filter to the specific intervals (EK and LK)
colls_EK <- colls %>% filter(epoch == "Early Cretaceous") 

# generate paleocoordinates (needs GPlates to be installed if using a Mac)
paleo <- reconstruct(colls_EK[, c("lng", "lat")], age=colls_EK$mid_ma, model="MERDITH2021", enumerate=FALSE)
colnames(paleo) <- c("plng", "plat")

# add to the rest
colls_EK<- cbind(colls_EK, paleo)

# subset species
hot_sp <- cloud_EK[which(cloud_EK$seasonal_temp>40),]
# assign EK age
hot_sp$age <- 121
# reconstruct model
palgeoEK_hot <- reconstruct("coastlines", age = 121, model="MERDITH2021")
## map theme
palaeomap_theme <- theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                                           axis.title.y=element_blank(), axis.text.y=element_blank(),
                                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                           legend.title=element_blank())
Map_hot_fx <-  ggplot() +
  ## Landmasses
  geom_sf(data = palgeoEK_hot, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = colls_EK, aes(x = plng, y = plat), color = "#FF8899", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("Hot species") +
  ## theme
  palaeomap_theme
Map_hot_fx


# seasonal temp LK
t_LK <- ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = seasonal_temp, fill = Pterosaur_taxa)) + 
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
t_LK

# Mann-Whitney U test
Wilcox_T <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_T_EK <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_T_LK <- wilcox.test(seasonal_temp ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## MAP EK
map_EK <- ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = MAP, fill = Pterosaur_taxa)) + 
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
map_EK

## MAP LK
map_LK <- ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = MAP, fill = Pterosaur_taxa)) + 
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
map_LK

# Mann-Whitney U test
Wilcox_MAP <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_MAP_EK <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_MAP_LK <- wilcox.test(MAP ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)

## seasonal Precip EK
p_EK <- ggplot(cloud_EK, aes(x = Pterosaur_taxa, y = seasonal_precip, fill = Pterosaur_taxa)) + 
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
p_EK
  
## seasonal Precip LK
p_LK <- ggplot(cloud_LK, aes(x = Pterosaur_taxa, y = seasonal_precip, fill = Pterosaur_taxa)) + 
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
p_LK

# Mann-Whitney U test
Wilcox_P <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_data, exact = FALSE)
Wilcox_P_EK <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_EK, exact = FALSE)
Wilcox_P_LK <- wilcox.test(seasonal_precip ~ Pterosaur_taxa, data = cloud_LK, exact = FALSE)



# plot as 2 grids
EK <- ggarrange(mat_EK, t_EK, map_EK, p_EK, ncol = 2, nrow = 2, 
                common.legend = TRUE, legend = "none")
LK <- ggarrange(mat_LK, t_LK, map_LK, p_LK, ncol = 2, nrow = 2, 
                common.legend = TRUE, legend = "none")
rain_all <- ggarrange(mat_EK, t_EK, map_EK, p_EK, 
                      mat_LK, t_LK, map_LK, p_LK, ncol = 4, nrow = 2, 
                      common.legend = TRUE, legend = "none") 

EK
LK
rain_all

ggsave(plot = rain_all,
       width = 14, height = 12, dpi = 600, 
       filename = "./plots/rain_all.pdf", useDingbats=FALSE)






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