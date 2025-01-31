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
  mutate_rows(group = PCA_data_earlyK$two.groups) ## FOR AZHD AND PTER $two.groups

summary(eK_pca)
eK_pca$rotation

confellip_eK <- eK_pca %>% 
  ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Early Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#FFA93D", "#572AAC" 
    #,"#248528","#D7E05A"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#FFA93D", "#572AAC" 
    #,"#248528","#D7E05A"
  )) ## add more colours if >n families!
confellip_eK

# Late Cretaceous
lK_pca <- PCA_data_lateK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_lateK$two.groups) ## FOR AZHD AND PTER $two.groups

summary(lK_pca)
lK_pca$rotation


confellip_lK <- lK_pca %>% 
  ordr::ggbiplot(data = PCA_data_lateK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Late Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#FFA93D", "#572AAC" 
  #  ,"#248528","#D7E05A","#993344","#0891A3","#1E4466"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#FFA93D", "#572AAC" 
  #  ,"#248528","#D7E05A", "#993344", "#0891A3","#1E4466"
  )) ## add more colours if >n families!

confellip_lK


# 2. PCA: statistics -------------------------------------------------------

# npMANOVA 
# https://www.rdocumentation.org/packages/RVAideMemoire/versions/0.9-80/topics/pairwise.perm.manova
# citation("RVAideMemoire")

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





########################################################################
# remove Lagerstatten effect by removing fossils from Lagerstatten
# Crato and Santana (Brazil), Solnhofen (Germany), Jehol (China)
#######################################################################





# Work in progress ------------------------------------------------------- 

####### Alternative for PCA
####################################### early Jurassic
## Change row names to occurrence numbers
PCA_data_earlyJ_df <- as.data.frame(PCA_data_earlyJ) # switch to dataframe
rownames(PCA_data_earlyJ_df) <- PCA_data_earlyJ_df[,1] # copy occurrence_no to row name
PCA_data_earlyJ_df$occurrence_no <- NULL # remove "occurrence_no" column
head(PCA_data_earlyJ_df) # check


## Pull out columns
earlyJ_pca <- PCA_data_earlyJ_df[, 1:4]
earlyJ_groups <- as.factor(PCA_data_earlyJ_df$family)


## PCA
earlyJ_pca_out <- prcomp(earlyJ_pca, center = TRUE, scale. = TRUE) 
summary(earlyJ_pca_out)

earlyJ_pca.var <- earlyJ_pca_out$sdev^2 # amount of variation each PC counts for
earlyJ_pca.var.per <- round(earlyJ_pca.var/sum(earlyJ_pca.var)*100, 1) # percentages

barplot(earlyJ_pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")


## get the name of the top 10 measurements that contribute most to pc1.
loading_scores <- earlyJ_pca_out$rotation[,1]
climate_scores <- abs(loading_scores) ## get the magnitudes
climate_score_ranked <- sort(climate_scores, decreasing=TRUE)
top_4 <- names(climate_score_ranked[1:4])

earlyJ_pca_out$rotation[top_4,1] ## show the scores (and +/- sign)

## Plot
earlyJ_PCA_plot <- ggbiplot::ggbiplot(earlyJ_pca_out, obs.scale = 1, var.scale = 1, groups = earlyJ_groups, ellipse = TRUE)
earlyJ_PCA_plot <- earlyJ_PCA_plot + scale_colour_manual(values = c("#0891A3", "#FFA93D", "#B00B69",  "#248528","#572AAC" )) + #"grey50",
  theme(panel.background = element_blank(),
        legend.position = "top", #legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA))
earlyJ_PCA_plot


###################################### middle Jurassic
## Change row names to occurrence numbers
PCA_data_midJ_df <- as.data.frame(PCA_data_midJ) # switch to dataframe
rownames(PCA_data_midJ_df) <- PCA_data_midJ_df[,1] # copy occurrence_no to row name
PCA_data_midJ_df$occurrence_no <- NULL # remove "occurrence_no" column
head(PCA_data_midJ_df) # check


## Pull out columns
midJ_pca <- PCA_data_midJ_df[, 1:4]
midJ_groups <- as.factor(PCA_data_midJ_df$family)


## PCA
midJ_pca_out <- prcomp(midJ_pca, center = TRUE, scale. = TRUE) 
summary(midJ_pca_out)

midJ_pca.var <- midJ_pca_out$sdev^2 # amount of variation each PC counts for
midJ_pca.var.per <- round(midJ_pca.var/sum(midJ_pca.var)*100, 1) # percentages

barplot(midJ_pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")


## get the name of the top 10 measurements that contribute most to pc1.
loading_scores <- midJ_pca_out$rotation[,1]
climate_scores <- abs(loading_scores) ## get the magnitudes
climate_score_ranked <- sort(climate_scores, decreasing=TRUE)
top_4 <- names(climate_score_ranked[1:4])

midJ_pca_out$rotation[top_4,1] ## show the scores (and +/- sign)

## Plot
midJ_PCA_plot <- ggbiplot::ggbiplot(midJ_pca_out, obs.scale = 1, var.scale = 1, groups = midJ_groups, ellipse = TRUE)
midJ_PCA_plot <- midJ_PCA_plot + scale_colour_manual(values = c("#0891A3", "#FFA93D", "#B00B69",  "#248528","#572AAC" )) + #"grey50",
  theme(panel.background = element_blank(),
        legend.position = "top", #legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA))
midJ_PCA_plot # error führende Minor der Ordnung 2 ist nicht positiv definit




# ------------------------------------------------------------------------------
### trying new groups: those that are phylogenetically useful
## Yu et al. 2023 tree 3 (most parsimonious) and (all) PBDB species

# for every "family" that is in the edited PBDB data (ptero_grouping_copy.csv)

### 1.
ptero_grouping[(which(ptero_grouping$family == "Campylognathoididae")),] 
# kick Campylognathoides_liasicus, Campylognathoides_zitteli
ptero_grouping[(which(ptero_grouping$family == "Pterosauria")),] 

# merge to one group, total of 6
A <- c("Austriadactylus_cristatus", "Caviramus_schesaplanensis", "Preondactylus_buffarinii",
       "Eudimorphodon_ranzii", "Eudimorphodon_rosenfeldi", "Raeticodactylus_filisurensis")

### 2.
ptero_grouping[(which(ptero_grouping$family == "Rhamphorhynchidae")),] 
# kick Sordes Pilosus, Orientognathus chaoyangensis,Pterorhynchus wellnhoferi

# one group, total of 8
B <- c("Angustinaripterus_longicephalus", "Cacibupteryx_caribensis", "Dorygnathus_banthensis",
       "Nesodactylus_hesperius", "Qinglongopterus_guoi", "Rhamphorhynchus_muensteri",
       "Scaphognathus_crassirostris", "Sericipterus_wucaiwanensis")

### 3.
ptero_grouping[(which(ptero_grouping$family == "Anurognathidae")),] 

# one group, total of 5
C <- c("Anurognathus_ammoni", "Batrachognathus_volans", "Dendrorhynchoides_curvidentatus",
       "Dendrorhynchoides_mutoudengensis", "Jeholopterus_ningchengensis")

### 4.
ptero_grouping[(which(ptero_grouping$family == "Pteranodontia")),] 
# put Pterodactylus_antiquus, Pterodactylus_kochi into "Gretheroidae"
ptero_grouping[(which(ptero_grouping$family == "Pterodactyloidea")),] 
# kick Kryptodrakon_progenitor, put all 3 Nyctosaurus into Pteranodontia, Ardeadactylus into "Gretheroidae"

# merge into one, total of 8
D <- c("Alamodactylus_byrdi", "Alcione_elainus", "Pteranodon_longiceps",
       "Pteranodon_sternbergi", "Simurghia_robusta", "Nyctosaurus_gracilis",
       "Nyctosaurus_lamegoi", "Nyctosaurus_nanus")

### 5.

ptero_grouping[(which(ptero_grouping$family == "Lophocratia")),] 
# kick Eosipterus_yangi
# + Protazhdarchidae (Aurorazhdarcho_primordius), + Macronychoptera (Cuspicephalus) + Liaodactylus (ctenochasmatidae)
# and Moganopterus (Boreopteridae)

# merge into one, total of 20
E <- c("Beipiaopterus_chenianus", "Ctenochasma_elegans", "Ctenochasma_porocristata",
       "Cycnorhamphus_suevicus", "Elanodactylus_prolatus", "Feilongus_youngi",
       "Gallodactylus_canjuersensis", "Germanodactylus_cristatus", "Germanodactylus_rhamphastinus",
       "Gnathosaurus_macrurus", "Gnathosaurus_subulatus", "Huanhepterus_quingyangensis",
       "Kepodactylus_insperatus", "Normannognathus_wellnhoferi", "Plataleorhynchus_streptophorodon",
       "Pterodaustro_guinazui", "Aurorazhdarcho_primordius", "Cuspicephalus_scarfi",
       "Liaodactylus_primus", "Moganopterus_zhuiana" )


# Macronychoptera
# Ctenochasmatidae
# Boreopteridae
# Ornithocheiroidea + 1 Borepteridae (Boreopterus_cuiae)
