# ******************************************************
#
#   Master Thesis
#
#   Mapping occurrences of pterosaurs
#
#   Carolin Grether
# ______________________________________________________
#
#   Modern map and Paleomaps
# 
# ******************************************************


# 0. Packages used in this script -----------------------------------------

library(tidyverse)
library(geoscale) # for plotting with the geological time scale on the x-axis (uses base R syntax)
library(viridis) # for colour scales
library(rgplates) # paleogeographic reconstructions
library(ggplot2) # for plotting
library(ggpubr) # for plotting
library(dplyr) # functions

select <- dplyr::select # ensure the select function is coming from dplyr

# ------------------------------------------------------------
# Data
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20)
occurrences_sp <- occurrences %>% filter(accepted_rank == "species")
## take necessary columns
taxon_inf <- select(occurrences_sp, collection_no, occurrence_no, accepted_name, 
                    early_interval, late_interval, max_ma, min_ma, formation)



ints_standard <- read.csv2("Data/Input/ints_standard_copy.csv") # manually fixed data
late_int <- read.csv2("Data/Input/stages_for_climate_data_late.csv")

# 0. put data into correct format -----------------------------------------------------------
## add stage names to data
Adat <- left_join(taxon_inf, ints_standard, by = "early_interval")
Bdat <- left_join(taxon_inf, late_int, by = "late_interval")

# merge datasets
big.dat <- rbind(Adat,Bdat)
# remove duplicates
taxon.dat <- big.dat[!duplicated(big.dat),]

# 1. Modern world map --------------------------------------------------------

## keep the info needed for making the map
locality_info <- occurrences_sp %>% 
  dplyr::select(collection_name, lat, lng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit()

## world map for ggplot to work with:
world_map <- map_data("world")
ggplot() + geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region)) 

## data
modern_map <- ggplot() + 
  geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region), 
           color = "grey80", fill = "grey90", linewidth = 0.1) +
  geom_point(data = locality_info, aes(lng, lat), alpha = 0.3, size = 4, colour = "#9B1999") +
  theme_void() + theme(legend.position = "none")
modern_map

## save as a pdf
ggsave(plot = modern_map,
       width = 8, height = 5, dpi = 600, 
       filename = "./plots/Modern_map.pdf", useDingbats=FALSE)



# 1. Palaeogeographic maps ------------------------------------------------

## simplified data object for palaeomap:
paleomap_data <- occurrences_sp %>% 
  select(collection_name, lat, lng, paleolat, paleolng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit(collection_name)

# add mean age to each collection
paleomap_data$mid_ma <- (paleomap_data$max_ma + paleomap_data$min_ma)/2


#-----------------------------------------------------------------------------
###### for mid ages (mean) ################
#-----------------------------------------------------------------------------
map_data_mid_MT <- paleomap_data %>% filter(mid_ma > 237.0) # no MT occurrence
map_data_mid_LT <- paleomap_data %>% filter(201.4 < mid_ma & mid_ma <= 237.0) 
map_data_mid_EJ <- paleomap_data %>% filter(174.7 < mid_ma & mid_ma <= 201.4)
map_data_mid_MJ <- paleomap_data %>% filter(161.5 < mid_ma & mid_ma <= 174.7)
map_data_mid_LJ <- paleomap_data %>% filter(145.0 < mid_ma & mid_ma <= 161.5)
map_data_mid_EK <- paleomap_data %>% filter(125.7 < mid_ma & mid_ma <= 145.0)
map_data_mid_LEK <- paleomap_data %>% filter(100.5 < mid_ma & mid_ma <= 125.7)
map_data_mid_ELK <- paleomap_data %>% filter(86.3 < mid_ma & mid_ma <= 100.5)
map_data_mid_LK <- paleomap_data %>% filter(66.0 < mid_ma & mid_ma <= 86.3)

## paleogeographies for the time bins from the GPlates (via rgplates)
paleogeog_LT <- reconstruct("static_polygons", age = 215, model="MERDITH2021") 
paleogeog_EJ <- reconstruct("coastlines", age = 190, model="MERDITH2021") 
paleogeog_MJ <- reconstruct("coastlines", age = 167, model="MERDITH2021")
paleogeog_LJ <- reconstruct("coastlines", age = 153, model="MERDITH2021")
paleogeog_EK <- reconstruct("coastlines", age = 135, model="MERDITH2021")
paleogeog_LEK <- reconstruct("coastlines", age = 113, model="MERDITH2021")
paleogeog_ELK <- reconstruct("coastlines", age = 93, model="MERDITH2021")
paleogeog_LK <- reconstruct("coastlines", age = 76, model="MERDITH2021")

## map theme
palaeomap_theme <- theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                                           axis.title.y=element_blank(), axis.text.y=element_blank(),
                                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                           legend.title=element_blank())


## Late Triassic
paleomap_mid_LT <-  ggplot() +
  ## andmasses
  geom_sf(data = paleogeog_LT, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_LT, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("Late Triassic") +
  ## theme
  palaeomap_theme
paleomap_mid_LT

## Early Jurassic
paleomap_mid_EJ <-  ggplot() +
  ## andmasses
  geom_sf(data = paleogeog_EJ, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_EJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("Early Jurassic") +
  ## theme
  palaeomap_theme
paleomap_mid_EJ

## Middle Jurassic
paleomap_mid_MJ <-  ggplot() +
  ## andmasses
  geom_sf(data = paleogeog_MJ, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_MJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ##  
  ggtitle("Middle Jurassic") +
  ## theme
  palaeomap_theme
paleomap_mid_MJ

## Late Jurassic
paleomap_mid_LJ <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LJ, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_LJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("Late Jurassic") +
  ## theme
  palaeomap_theme
paleomap_mid_LJ

## Early Cretaceous
paleomap_mid_EK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_EK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_EK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("early Early Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_mid_EK

## Late Early Cretaceous
paleomap_mid_LEK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LEK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_LEK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("late Early Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_mid_LEK

## Early Late Cretaceous
paleomap_mid_ELK <-  ggplot() +
  ## andmasses
  geom_sf(data = paleogeog_ELK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_ELK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("early Late Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_mid_ELK

## Late Cretaceous
paleomap_mid_LK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_mid_LK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("late Late Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_mid_LK

## save as a .pdf
ggsave(plot = paleomap_mid_LT,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_LateTriassic.pdf", useDingbats=FALSE)
# EJ
ggsave(plot = paleomap_mid_EJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_EarlyJurassic.pdf", useDingbats=FALSE)
# MJ
ggsave(plot = paleomap_mid_MJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_MiddleJurassic.pdf", useDingbats=FALSE)
# LJ
ggsave(plot = paleomap_mid_LJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_LateJurassic.pdf", useDingbats=FALSE)
# EK
ggsave(plot = paleomap_mid_EK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_erEarlyCretaceous.pdf", useDingbats=FALSE)
# LEK
ggsave(plot = paleomap_mid_LEK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_ltEarlyCretaceous.pdf", useDingbats=FALSE)
# ELK
ggsave(plot = paleomap_mid_ELK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_erLateCretaceous.pdf", useDingbats=FALSE)
# LK
ggsave(plot = paleomap_mid_LK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_mid_ltLateCretaceous.pdf", useDingbats=FALSE)

# plot as 2 grids
T_and_J_mid <- ggarrange(paleomap_mid_LT, paleomap_mid_EJ, paleomap_mid_MJ, paleomap_mid_LJ, ncol = 2, nrow = 2)
K_in_4_mid <- ggarrange(paleomap_mid_EK, paleomap_mid_LEK, paleomap_mid_ELK, paleomap_mid_LK, ncol = 2, nrow = 2)
Mid_all <- ggarrange(paleomap_mid_LT, paleomap_mid_EJ, paleomap_mid_MJ, paleomap_mid_LJ, 
                     paleomap_mid_EK, paleomap_mid_LEK, paleomap_mid_ELK, paleomap_mid_LK, ncol = 2, nrow = 4) 

T_and_J_mid
K_in_4_mid
Mid_all

ggsave(plot = Mid_all,
       width = 11, height = 14, dpi = 600, 
       filename = "./plots/Palaeomap_mid_all.pdf", useDingbats=FALSE)




#-----------------------------------------------------------------------------
###### only FAD (earliest occurences) ################
#-----------------------------------------------------------------------------


map_data_MT <- paleomap_data %>% filter(max_ma > 237.0) # one MT occurrence
map_data_LT <- paleomap_data %>% filter(201.4 < max_ma & max_ma <= 237.0) 
map_data_EJ <- paleomap_data %>% filter(174.7 < max_ma & max_ma <= 201.4)
map_data_MJ <- paleomap_data %>% filter(161.5 < max_ma & max_ma <= 174.7)
map_data_LJ <- paleomap_data %>% filter(145.0 < max_ma & max_ma <= 161.5)
map_data_EK <- paleomap_data %>% filter(125.7 < max_ma & max_ma <= 145.0)
map_data_LEK <- paleomap_data %>% filter(100.5 < max_ma & max_ma <= 125.7)
map_data_ELK <- paleomap_data %>% filter(86.3 < max_ma & max_ma <= 100.5)
map_data_LK <- paleomap_data %>% filter(66.0 <= max_ma & max_ma <= 86.3)


## Late Triassic
paleomap_LT <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LT, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_LT, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## add MT occurrence
  geom_point(data = map_data_MT, aes(x = paleolng, y = paleolat), color = "#FFA93D", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("Late Triassic") +
  ## theme
  palaeomap_theme
paleomap_LT


## Early Jurassic
paleomap_EJ <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_EJ, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_EJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("Early Jurassic") +
  ## theme
  palaeomap_theme
paleomap_EJ

## Middle Jurassic
paleomap_MJ <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_MJ, colour = "grey75", fill = "grey75") +
  ## occurrence data 
  geom_point(data = map_data_MJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("Middle Jurassic") +
  ## theme
  palaeomap_theme
paleomap_MJ

## Late Jurassic
paleomap_LJ <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LJ, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_LJ, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("Late Jurassic") +
  ## theme
  palaeomap_theme
paleomap_LJ

## early Early Cretaceous
paleomap_EK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_EK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_EK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("early Early Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_EK

## late Early Cretaceous
paleomap_LEK <-  ggplot() +
  ## andmasses
  geom_sf(data = paleogeog_LEK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_LEK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("late Early Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_LEK

## early Late Cretaceous
paleomap_ELK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_ELK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_ELK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title
  ggtitle("early Late Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_ELK

## late Late Cretaceous
paleomap_LK <-  ggplot() +
  ## landmasses
  geom_sf(data = paleogeog_LK, colour = "grey75", fill = "grey75") +
  ## occurrence data
  geom_point(data = map_data_LK, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
  ## title 
  ggtitle("late Late Cretaceous") +
  ## theme
  palaeomap_theme
paleomap_LK


## save as a .pdf
# LT
ggsave(plot = paleomap_LT,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_LateTriassic.pdf", useDingbats=FALSE)
# EJ
ggsave(plot = paleomap_EJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_EarlyJurassic.pdf", useDingbats=FALSE)
# MJ
ggsave(plot = paleomap_MJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_MiddleJurassic.pdf", useDingbats=FALSE)
# LJ
ggsave(plot = paleomap_LJ,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_LateJurassic.pdf", useDingbats=FALSE)
# EK
ggsave(plot = paleomap_EK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_EarlyCretaceous.pdf", useDingbats=FALSE)
# LEK
ggsave(plot = paleomap_LEK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_ltEarlyCretaceous.pdf", useDingbats=FALSE)
# ELK
ggsave(plot = paleomap_ELK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_erLateCretaceous.pdf", useDingbats=FALSE)
# LK
ggsave(plot = paleomap_LK,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomap_LateCretaceous.pdf", useDingbats=FALSE)

# plot as 2 grids
T_and_J <- ggarrange(paleomap_LT, paleomap_EJ, paleomap_MJ, paleomap_LJ, ncol = 2, nrow = 2)
K_in_4 <- ggarrange(paleomap_EK, paleomap_LEK, paleomap_ELK, paleomap_LK, ncol = 2, nrow = 2)
FAD_all <- ggarrange(paleomap_LT, paleomap_EJ, paleomap_MJ, paleomap_LJ, 
          paleomap_EK, paleomap_LEK, paleomap_ELK, paleomap_LK, ncol = 2, nrow = 4) 

T_and_J
K_in_4
FAD_all

ggsave(plot = FAD_all,
       width = 11, height = 14, dpi = 600, 
       filename = "./plots/Palaeomap_FAD_all.pdf", useDingbats=FALSE)
