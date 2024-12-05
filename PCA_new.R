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

library(devtools) # to install ggord
install_github("fawda123/ggord") # For PCA plots
library(ggord)
library(ordr) # manipulating data objects
library(ggplot2)
library(ggdist)


# 1. ALL DATA (Yu et al.2023) Organise data -----------------------------------------------------------

## Import species/climate data if not already loaded:
species_climate <- read.csv("Data/climate/species_climate.csv")
glimpse(species_climate)

# groups for PCA
ptero_grouping <- read.csv2("Data/Input/ptero_groups_copy.csv") # all PBDB pterosaurs
ptero_grouping <- read.csv2("Data/Input/azhd_and_pteran.csv") # ONLY AZHDARCHOIDEA AND PTERANODONTIA


## merge species climate and group data here
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

unique(species_climate_group$family)

## Only the columns needed for the PCA:
# PCA_data <- subset(species_climate_group, select = c(occurrence_no,  ## ALL PTEROS
#                                                      MAT, seasonal_temp,
#                                                      MAP, seasonal_precip,
#                                                      family, ## ***** THIS IS THE GROUPING VARIABLE ****
#                                                      early_interval
# ))

PCA_data <- subset(species_climate_group, select = c(occurrence_no,   ## ONLY FOR AZHD AND PTERAN
                                                     MAT, seasonal_temp,
                                                     MAP, seasonal_precip,
                                                     two.groups, ## ***** THIS IS THE GROUPING VARIABLE ****
                                                     early_interval
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

## Filter to the 4 specific intervals
## For Azhd and Pteran, Cretaceous:
PCA_data_earlyK <- PCA_data %>% filter(epoch == "Early Cretaceous") 
PCA_data_lateK <- PCA_data %>% filter(epoch == "Late Cretaceous") 

## REST (FOR ALL PTEROS)
PCA_data_earlyJ <- PCA_data %>% filter(epoch == "Early Jurassic") 
PCA_data_midJ <- PCA_data %>% filter(epoch == "Middle Jurassic") 
PCA_data_lateJ <- PCA_data %>% filter(epoch == "Late Jurassic") 

unique(PCA_data_all_lateK$family) # check all groups, best is n<6


####### ALL PTEROS: merge some groups in Late K

 PCA_all_lateK_prun <- PCA_data_lateK[-c(which(PCA_data_lateK$family == "Lonchodraconidae")),]
 PCA_all_lateK_prun <- PCA_all_lateK_prun[-c(which(PCA_all_lateK_prun$family == "Nyctosauridae")),]
 PCA_all_lateK_prun <- PCA_all_lateK_prun[-c(which(PCA_all_lateK_prun$family == "Azhdarchoidea")),]

# keep
 PCA_concise_LK <- PCA_all_lateK_prun
 
 # Pteranodontoidea
 PCA_concise_LK$family[which(PCA_concise_LK$family=="Pteranodontia")] <- "Pteranodontoidea"
 PCA_concise_LK$family[which(PCA_concise_LK$family=="Pteranodontidae")] <- "Pteranodontoidea"
 PCA_con_LK <- PCA_concise_LK
 
 # Ornithocheiroidea
 PCA_concise_LK$family[which(PCA_concise_LK$family=="Ornithocheiridae")] <- "Ornithocheiroidea"



# PCA: analysis and plots --------------------------------------------------

## Perform a separate PCA and plot for each interval


# early Jurassic
eJ_pca <- PCA_data_earlyJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_earlyJ$family)

summary(eJ_pca)
eJ_pca$rotation

confellip_eJ <- eJ_pca %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) + ## add more colours if >n families!
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) ## add more colours if >n families!
confellip_eJ # error führende Minor der Ordnung 2 ist nicht positiv definit


# middle Jurassic
mJ_pca <- PCA_data_midJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_midJ$family)

summary(mJ_pca)
mJ_pca$rotation 

confellip_mJ <- mJ_pca %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) + ## add more colours if >n families!
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) ## add more colours if >n families!
confellip_mJ


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

# LK concise Ornithocheiridae + Ornithocheiroidea
lK_pca_con <- PCA_con_LK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_con_LK$family) ## FOR AZHD AND PTER $two.groups

confellip_lK_con <- lK_pca_con %>% 
  ordr::ggbiplot(data = PCA_con_LK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Late Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    ,"#248528","#D7E05A","#993344"
    #"#0891A3","#1E4466"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    ,"#248528","#D7E05A", "#993344"
    #"#0891A3","#1E4466"
  )) ## add more colours if >n families!

confellip_lK_con

# LK concise, only Ornithocheiroidea
lK_pca_concise <- PCA_concise_LK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_concise_LK$family) ## FOR AZHD AND PTER $two.groups

confellip_lK_concise <- lK_pca_concise %>% 
  ordr::ggbiplot(data = PCA_concise_LK ,aes(color = group)) +
  theme_bw() +
  ggtitle("Late Cretaceous") +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    ,"#248528","#D7E05A"
    #,"#993344"
    #"#0891A3","#1E4466"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    ,"#248528","#D7E05A"
    #, "#993344"
    #"#0891A3","#1E4466"
  )) ## add more colours if >n families!

confellip_lK_concise


# UNCOMMENT FOLLOWING FOR VENDITTI DATA
# Organise data -----------------------------------------------------------------------------------

## Import species/climate data if not already loaded:
species_climate <- read.csv("Data/climate/species_climate.csv")

## Venditti data
Ven.dat <- read.csv2("Data/Output/Species_eff_data.csv")
# To include Quetzalcoatlus data, change eff.dat "Quetzalcoatlus spp" to "Quetzalcoatlus northropi"
Ven.dat$species[62] <- "Quetzalcoatlus_northropi"
Ven.dat$species[36] <- "Hatzegopteryx_thambema" # fixed spelling
Ven.dat$species[40] <- "Huaxiapterus_corollatus" # fixed spelling

# groups for PCA
ptero_grouping <- read.csv2("Data/Input/ptero_groups_copy.csv")

# merge groupings and climate
Ven_taxa <- merge(ptero_grouping, Ven.dat, by.x = "ptero_taxa", by.y = "species")
Ven_climate <- merge(species_climate, Ven_taxa, by.x = "accepted_name", by.y = "ptero_taxa")


## Only the columns needed for the PCA:
PCA_Ven <- subset(Ven_climate, select = c(occurrence_no,  ## ALL PTEROS
                                          MAT, seasonal_temp, 
                                          MAP, seasonal_precip, 
                                          family, ## ***** THIS IS THE GROUPING VARIABLE **** 
                                          early_interval
))

## Load auxillary dataset to standardise the time intervals
ints_standard <- read_csv2("Data/Input/ints_standard_copy.csv") # manually fixed data

## rename 'early interval'
PCA_Ven <- rename(PCA_Ven, interval_std = early_interval)
ints_standard <- rename(ints_standard, interval_std = early_interval)

## join aux. dataset to main PCA data
PCA_Ven <- left_join(PCA_Ven, ints_standard, by = "interval_std")

## Remove those that have NA (stratigraphic range is too long)
PCA_Ven <- na.omit(PCA_Ven, epoch)

## Filter to the specific intervals
PCA_Ven_eJ <- PCA_Ven %>% filter(epoch == "Early Jurassic") 
PCA_Ven_mJ <- PCA_Ven %>% filter(epoch == "Middle Jurassic") 
PCA_Ven_lJ <- PCA_Ven %>% filter(epoch == "Late Jurassic") 

# Cretaceous in 2
PCA_Ven_eK <- PCA_Ven %>% filter(epoch == "Early Cretaceous") 
PCA_Ven_lK <- PCA_Ven %>% filter(epoch == "Late Cretaceous") 

# Cretaceous in 4
K_one <- c("Berriasian", "Valanginian", "Hauterivian", "Barremian")
K_two <- c("Aptian", "Albian")
K_three <- c("Cenomanian", "Turonian", "Coniacian", "Santonian")
K_four <- c("Campanian", "Maastrichtian")

PCA_Ven_K1 <- PCA_Ven[PCA_Ven$stage%in%K_one,] # CHANGE CODE to PCA_data if checking "all"
PCA_Ven_K2 <- PCA_Ven[PCA_Ven$stage%in%K_two,] 
PCA_Ven_K3 <- PCA_Ven[PCA_Ven$stage%in%K_three,]
PCA_Ven_K4 <- PCA_Ven[PCA_Ven$stage%in%K_four,]

unique(PCA_Ven_K1$family) # check all groups, best is n<6

## PCA
# early Jurassic
Ven_eJ_pca <- PCA_Ven_eJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_eJ$family)

summary(Ven_eJ_pca)
Ven_eJ_pca$rotation

confellip_Ven_eJ <- Ven_eJ_pca %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528")) + ## add more colours if >n families! #, "#FFA93D","#572AAC"
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528")) ## add more colours if >n families!
confellip_Ven_eJ # error führende Minor der Ordnung 2 ist nicht positiv definit


# middle Jurassic
Ven_mJ_pca <- PCA_Ven_mJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_mJ$family)

summary(Ven_mJ_pca)
Ven_mJ_pca$rotation 

confellip_Ven_mJ <- Ven_mJ_pca %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(  "#248528", "#FFA93D","#572AAC")) + ## add more colours if >n families! #"#0891A3","#1E44AA",
  scale_fill_manual(values = c(  "#248528", "#FFA93D","#572AAC")) ## add more colours if >n families!
confellip_Ven_mJ


# late Jurassic
Ven_lJ_pca <- PCA_Ven_lJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_lJ$family)

summary(Ven_lJ_pca)
Ven_lJ_pca$rotation 

confellip_Ven_lJ <- Ven_lJ_pca %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) + ## add more colours if >n families! #"#0891A3","#1E44AA",
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) ## add more colours if >n families!
confellip_Ven_lJ


# early Cretaceous
eK_pca <- PCA_data_eK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_earlyK$family) ## FOR AZHD AND PTER $two.groups

summary(eK_pca)
eK_pca$rotation

confellip_eK <- eK_pca %>% 
  ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    #,"#248528","#D7E05A"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    #,"#248528","#D7E05A"
  )) ## add more colours if >n families!
confellip_eK

# Late Cretaceous
lK_pca <- PCA_data_lateK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_lateK$family) ## FOR AZHD AND PTER $two.groups

confellip_lK <- lK_pca %>% 
  ordr::ggbiplot(data = PCA_data_lateK ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    #,"#248528","#D7E05A"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D"
    #,"#248528","#D7E05A"
  )) ## add more colours if >n families!

confellip_lK

################ Cretaceous in 4 slices

# K 1
K1_pca <- PCA_Ven_K1[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_K1$family) ## FOR AZHD AND PTER $two.groups

summary(K1_pca)
K1_pca$rotation

confellip_K1 <- K1_pca %>% 
  ordr::ggbiplot(data = PCA_Ven_K1 ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D","#248528"
    ,"#248528","#D7E05A", "#FFFF00", "#993344"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA",
    "#572AAC" , "#FFA93D","#248528"
    ,"#248528","#D7E05A","#FFFF00", "#993344"
  )) ## add more colours if >n families!
confellip_K1

# K 2
K2_pca <- PCA_Ven_K2[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_K2$family) ## FOR AZHD AND PTER $two.groups

summary(K2_pca)
K2_pca$rotation

confellip_K2 <- K2_pca %>% 
  ordr::ggbiplot(data = PCA_Ven_K2 ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA",
                                 "#572AAC" , "#FFA93D","#248528"
                                 ,"#248528","#D7E05A", "#FFFF00","#FF1199", "#993344"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA",
                               "#572AAC" , "#FFA93D","#248528"
                               ,"#248528","#D7E05A", "#FFFF00","#FF1199", "#993344"
  )) ## add more colours if >n families!
confellip_K2

# K3
K3_pca <- PCA_Ven_K3[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_K3$family) ## FOR AZHD AND PTER $two.groups

summary(K3_pca)
K3_pca$rotation

confellip_K3 <- K3_pca %>% 
  ordr::ggbiplot(data = PCA_Ven_K3 ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA",
                                 "#572AAC" , "#FFA93D","#248528"
                                 ,"#248528","#D7E05A", "#FFFF00","#993344"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA",
                               "#572AAC" , "#FFA93D","#248528"
                               ,"#248528","#D7E05A", "#FFFF00","#993344"
  )) ## add more colours if >n families!
confellip_K3

# K4
K4_pca <- PCA_Ven_K4[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_Ven_K4$family) ## FOR AZHD AND PTER $two.groups

summary(K4_pca)
K4_pca$rotation

confellip_K4 <- K4_pca %>% 
  ordr::ggbiplot(data = PCA_Ven_K4 ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c("#0891A3", "#1E44AA",
                                 "#572AAC" , "#FFA93D","#248528"
                                 #,"#248528","#D7E05A", "#FFFF00"
  )) + ## add more colours if >n families!
  scale_fill_manual(values = c("#0891A3", "#1E44AA",
                               "#572AAC" , "#FFA93D","#248528"
                               #,"#248528","#D7E05A", "#FFFF00"
  )) ## add more colours if >n families!
confellip_K4


########################################################################
# remove Lagerstatten effect by removing fossils from Lagerstatten
# Crato and Santana (Brazil), Solnhofen (Germany), Jehol (China)
#######################################################################




###############################################################################
###############################################################################
#-----------------------------------------------------------------------------
# RAINCLOUD PLOTS
# for Pteranodontia and Azhdarchoidea

#-----------------------------------------------------------------------------
###############################################################################
###############################################################################

# rename data
cloud_data <- PCA_data
colnames(cloud_data)[6] <- "Pterosaur_taxa"
colnames(cloud_data)

# one dataset each for Azhd and Pteran
# cloud_az <- cloud_data %>% 
#   arrange(Pterosaur_taxa) %>% 
#   .[which(.$Pterosaur_taxa=="Azhdarchoidea"),]
# 
# cloud_pt <- cloud_data %>% 
#   arrange(Pterosaur_taxa) %>% 
#   .[which(.$Pterosaur_taxa=="Pteranodontoidea"),]


# Raincloud plots - entire Cretaceous
## MAT
ggplot(cloud_data, aes(x = Pterosaur_taxa, y = MAT, fill = Pterosaur_taxa)) + 
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


## seasonal temp
ggplot(cloud_data, aes(x = Pterosaur_taxa, y = seasonal_temp, fill = Pterosaur_taxa)) + 
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


## MAP
ggplot(cloud_data, aes(x = Pterosaur_taxa, y = MAP, fill = Pterosaur_taxa)) + 
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
  labs(x = NULL, y = "Mean Annual Precipitation") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")

## seasonal Precip
ggplot(cloud_data, aes(x = Pterosaur_taxa, y = seasonal_precip, fill = Pterosaur_taxa)) + 
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
  labs(x = NULL, y = "Seasonal Precipitation") + #coord_flip()
  coord_cartesian(xlim = c(1.2, NA), clip = "off")


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
