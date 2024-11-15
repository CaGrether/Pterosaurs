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


# 1. Organise data -----------------------------------------------------------

## Import species/climate data if not already loaded:
species_climate <- read.csv("Data/climate/species_climate.csv")
glimpse(species_climate)

# groups for PCA
ptero_grouping <- read.csv2("Data/Input/ptero_groups_copy.csv") # JURASSIC AND CRETACEOUS
#ptero_grouping <- read.csv2("Data/Input/azhd_and_pteran.csv") # ONLY AZHDARCHOIDEA AND PTERANODONTIA


## merge species climate and group data here
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

unique(species_climate_group$family)

## Only the columns needed for the PCA:
PCA_data <- subset(species_climate_group, select = c(occurrence_no,  ## ALL PTEROS
                                                     MAT, seasonal_temp, 
                                                     MAP, seasonal_precip, 
                                                     family, ## ***** THIS IS THE GROUPING VARIABLE **** 
                                                     early_interval
))

# PCA_data <- subset(species_climate_group, select = c(occurrence_no,   ## ONLY FOR AZHD AND PTERAN
#                                                      MAT, seasonal_temp, 
#                                                      MAP, seasonal_precip, 
#                                                      two.groups, ## ***** THIS IS THE GROUPING VARIABLE **** 
#                                                      early_interval
# ))


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

unique(PCA_data_lateK$family) # check all groups, best is n<6


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


# UNCOMMENT FOLLOWING FOR VENDITTI DATA
# ------------------------------------------------------------------------

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
PCA_Ven_eK <- PCA_Ven %>% filter(epoch == "Early Cretaceous") 
PCA_Ven_lK <- PCA_Ven %>% filter(epoch == "Late Cretaceous") 

unique(PCA_Ven_lJ$family) # check all groups, best is n<6

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
