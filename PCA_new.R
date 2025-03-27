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
#   PCA plots
# 
# ******************************************************


## Load package(s):
library(tidyverse)
library(FactoMineR) # For PCA plots

library(devtools) # to install ggord
install_github("fawda123/ggord") # For PCA plots
library(ggord)
library(ordr) # manipulating data objects
library(ggplot2)
library(ggdist)
library(RVAideMemoire) # for PCA statistics


# 0. Data preparation

## Import species/climate data if not already loaded:
species_climate <- read.csv("Data/climate/species_climate.csv")
glimpse(species_climate)

# groups for PCA
ptero_grouping <- read.csv2("Data/Input/azhd_and_pteran.csv") # Azhdarchoidea and Pteranodontoidea


## merge species climate and group data here
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

unique(species_climate_group$family)


PCA_data <- subset(species_climate_group, select = c(occurrence_no,   
                                                     MAT, seasonal_temp,
                                                     MAP, seasonal_precip,
                                                     two.groups, 
                                                     early_interval,
                                                     accepted_name
))


## Load auxillary dataset to standardise the time intervals
ints_standard <- read_csv2("Data/Input/ints_standard_copy.csv") # manually fixed data

## rename 'early interval'
PCA_data <- rename(PCA_data, interval_std = early_interval)
ints_standard <- rename(ints_standard, interval_std = early_interval)

## join aux. dataset to main PCA data
PCA_data <- left_join(PCA_data, ints_standard, by = "interval_std")

## Remove those that have NA (stratigraphic range is too long)
PCA_data <- na.omit(PCA_data, epoch)

## Filter to the specific intervals
## Cretaceous:
PCA_data_earlyK <- PCA_data %>% filter(epoch == "Early Cretaceous") 
PCA_data_lateK <- PCA_data %>% filter(epoch == "Late Cretaceous") 



# 1. PCA: analysis and plots --------------------------------------------------

# early Cretaceous
eK_pca <- PCA_data_earlyK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_earlyK$two.groups) 

summary(eK_pca)
eK_pca$rotation

confellip_eK <- eK_pca %>% 
  ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Early Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(
    "#FFA93D", "#572AAC" 
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(
    "#FFA93D", "#572AAC" 
  )) ## add more colours if >n families!
confellip_eK

# Late Cretaceous
lK_pca <- PCA_data_lateK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_lateK$two.groups)

summary(lK_pca)
lK_pca$rotation


confellip_lK <- lK_pca %>% 
  ordr::ggbiplot(data = PCA_data_lateK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Late Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(
    "#FFA93D", "#572AAC" 
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(
    "#FFA93D", "#572AAC" 
  )) ## add more colours if >n families!

confellip_lK


# 2. PCA: statistics -------------------------------------------------------

# The "fact" argument will be the groups (here: Azhdarchoidea and Pteranodontoidea)

EK_PCgroups <- PCA_data_earlyK$two.groups
head(EK_PCgroups)


LK_PCgroups <- PCA_data_lateK$two.groups
head(LK_PCgroups)

# The "resp" argument will be the euclidian distance of the PC scores. 

EK_PCscores <- as.data.frame(eK_pca$x) # PC scores from above
EK_PCscores <- cbind(EK_PCgroups, EK_PCscores)

LK_PCscores <- as.data.frame(lK_pca$x) # PC scores from above
LK_PCscores <- cbind(LK_PCgroups, LK_PCscores)


## npMANOVA:
npmanova_EK <- pairwise.perm.manova(dist(EK_PCscores[,1:2], "euclidian"), EK_PCscores$EK_PCgroups, nperm = 10000, p.method = "BH")
npmanova_LK <- pairwise.perm.manova(dist(LK_PCscores[,1:2], "euclidian"), LK_PCscores$LK_PCgroups, nperm = 10000, p.method = "BH")


# If p is lower than 0.05, the groups are significantly different (i.e., the PC scores of them are significantly different)

