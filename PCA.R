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
species_climate <- read.csv("Data/climate/species_climate.csv")
glimpse(species_climate)

ptero_grouping <- read.csv2("Data/Input/ptero_groups_copy.csv") # UNCOMMENT TO GET ALL PTEROS
#ptero_grouping <- read.csv2("Data/Input/azhd_and_pteran.csv") # ONLY AZHDARCHOIDEA AND PTERANODONTIA

## merge species climate and group data here
# only take species that we're using
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)
species_climate_group <- merge(species_climate, ptero_grouping, by.x = "accepted_name", by.y = "ptero_taxa")

unique(species_climate_group$family)

## Take only the columns we need for the PCA:
 PCA_data <- subset(species_climate_group, select = c(occurrence_no,  ## UNCOMMENT FOR ALL PTEROS
                                                MAT, seasonal_temp, 
                                                MAP, seasonal_precip, 
                                                family, ## ***** THIS WILL BE YOUR GROUPING VARIABLE **** 
                                                early_interval
 ))

#PCA_data <- subset(species_climate_group, select = c(occurrence_no,   ## ONLY FOR AZHD AND PTERAN
#                                                     MAT, seasonal_temp, 
#                                                     MAP, seasonal_precip, 
#                                                     # family,
#                                                     two.groups, ## ***** THIS WILL BE YOUR GROUPING VARIABLE **** 
#                                                     early_interval
#))


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

## Filter to the 4 specific intervals
## For Azhd and Pteran, Cretaceous:
PCA_data_earlyK <- PCA_data %>% filter(epoch == "Early Cretaceous") 
PCA_data_lateK <- PCA_data %>% filter(epoch == "Late Cretaceous") 

## UNCOMMENT FOR ALL PTEROS
PCA_data_earlyJ <- PCA_data %>% filter(epoch == "Early Jurassic") 
PCA_data_midJ <- PCA_data %>% filter(epoch == "Middle Jurassic") 
PCA_data_lateT <- PCA_data %>% filter(epoch == "Late Triassic") 
PCA_data_lateJ <- PCA_data %>% filter(epoch == "Late Jurassic") 

unique(PCA_data_earlyK$family) # check all groups, best is n<6

# PCA: analysis and plots --------------------------------------------------

## Perform a separate PCA and plot for each interval and 

# late Triassic
lT_pca <- PCA_data_lateT[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
mutate_rows(group = PCA_data_lateT$family) ## 

summary(lT_pca) # 1 explains ~73%, 2 explains ~17%
lT_pca$rotation
# in late Triassic: the variance of the (environmental?) data is explained best by precipitation
# while temperature is less important?

confellip_lT <- lT_pca %>% #na.omit(diet_group) %>%
  ordr::ggbiplot(aes(color = group)) +
  #theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#1E44AA", "#248528", "#FFA93D")) + ## add more colours if >5 families!
  scale_fill_manual(values = c( "#1E44AA", "#248528", "#FFA93D")) ## add more colours if >5 families!
confellip_lT

### for ProgPal poster
# confellip_lT <- lT_pca %>% #na.omit(diet_group) %>% 
#   ordr::ggbiplot(aes(color = group)) +
#   #theme_bw() +
#   theme(plot.background = element_rect(fill = '#ECDDBF'))+
#   ordr::geom_rows_point() +
#   geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
#   ordr::geom_cols_vector(color = "#444444") + # adds the arrows
#   scale_colour_manual(values = c( "#63c77cff", "#1E44AA", "#f5950fff")) + ## add more colours if >5 families!
#   scale_fill_manual(values = c( "#63c77cff", "#1E44AA", "#f5950fff")) ## add more colours if >5 families!
# confellip_lT


# early Jurassic
eJ_pca <- PCA_data_earlyJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_earlyJ$family)

summary(eJ_pca)
eJ_pca$rotation

confellip_eJ <- eJ_pca %>% #na.omit(diet_group) %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) + ## add more colours if >5 families!
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) ## add more colours if >5 families!
confellip_eJ


# middle Jurassic
mJ_pca <- PCA_data_midJ[,2:5] %>%
  prcomp(scale = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
  mutate_rows(group = PCA_data_midJ$family)

summary(mJ_pca)
mJ_pca$rotation 
# does this make any sense? seasonal T and MAP have high influence, MAT and seasonal P not?

confellip_mJ <- mJ_pca %>% #na.omit(diet_group) %>% 
  ordr::ggbiplot(aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows
  scale_colour_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) + ## add more colours if >5 families!
  scale_fill_manual(values = c( "#0891A3","#1E44AA", "#248528", "#FFA93D","#572AAC")) ## add more colours if >5 families!
confellip_mJ


# early Cretaceous
eK_pca <- PCA_data_earlyK[,2:5] %>%
  prcomp(center = T, scale. = TRUE) %>%
  ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
mutate_rows(group = PCA_data_earlyK$two.groups) ## FOR ALL PTEROS $family

summary(eK_pca)
eK_pca$rotation # seasonal temp

## Code for getting loadings
loadings_scores_1 <- eK_pca$rotation[,1]  # write 2 for PC2
climate_scores <- abs(loadings_scores_1) ## get the magnitudes
climate_score_ranked <- sort(climate_scores, decreasing=TRUE)
eK_pca$rotation[climate_score_ranked[1]] ## show the scores (and +/- sign)
##### THIS MAKES NO SENSE CHECK OVER
which(eK_pca$rotation==climate_score_ranked[1])


confellip_eK <- eK_pca %>% #na.omit(diet_group) %>%
  ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = group)) +
  theme_bw() +
  ordr::geom_rows_point() +
  geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
  ordr::geom_cols_vector(color = "#444444") + # adds the arrows, but where are the labels?
  scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
                                 "#572AAC" , "#FFA93D"
                                 #,"#248528","#D7E05A"
                                 )) + ## add more colours if >5 families!
  scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
                               "#572AAC" , "#FFA93D"
                               #,"#248528","#D7E05A"
                               )) ## add more colours if >5 families!

#### for ProgPal poster
# confellip_eK <- eK_pca %>% #na.omit(diet_group) %>% 
#   ordr::ggbiplot(data = PCA_data_earlyK ,aes(color = group)) +
#   #theme_bw() +
#   theme(plot.background = element_rect(fill = '#ECDDBF'))+
#   ordr::geom_rows_point() +
#   geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
#   ordr::geom_cols_vector(color = "#444444") + # adds the arrows
#   scale_colour_manual(values = c("#0891A3", "#63c77cff","#572AAC" , "#1E44AA",
#                                  "#f5950fff","#D7E05A")) + ## add more colours if >5 families!
#   scale_fill_manual(values = c("#0891A3", "#63c77cff","#572AAC" , "#1E44AA",
#                                "#f5950fff","#D7E05A")) ## add more colours if >5 families!
 confellip_eK
 
 # Late Cretaceous
 lK_pca <- PCA_data_lateK[,2:5] %>%
   prcomp(center = T, scale. = TRUE) %>%
   ordr::as_tbl_ord() %>% # if error, check package 'ordr' is installed correctly
   mutate_rows(group = PCA_data_lateK$two.groups) ## FOR ALL PTEROS $family
 
 confellip_lK <- lK_pca %>% #na.omit(diet_group) %>%
   ordr::ggbiplot(data = PCA_data_lateK ,aes(color = group)) +
   theme_bw() +
   ordr::geom_rows_point() +
   geom_polygon(aes(fill = group), color = NA, alpha = .25, stat = "rows_ellipse") +
   ordr::geom_cols_vector(color = "#444444") + # adds the arrows, but where are the labels?
   scale_colour_manual(values = c(#"#0891A3", "#1E44AA",
      "#572AAC" , "#FFA93D"
     #,"#248528","#D7E05A"
   )) + ## add more colours if >5 families!
   scale_fill_manual(values = c(#"#0891A3", "#1E44AA",
     "#572AAC" , "#FFA93D"
     #,"#248528","#D7E05A"
   )) ## add more colours if >5 families!
 
confellip_lK
 
 ## get new groups for PCA
 groupings_help <- subset(species_climate_group, select = c(occurrence_no,
                                                              family, ## ***** THIS WILL BE YOUR GROUPING VARIABLE **** 
                                                              early_interval,
                                                              accepted_name
 ))
 
 groupings_help <- rename(groupings_help, interval_std = early_interval)
 groupings_help <- left_join(groupings_help, ints_standard, by = "interval_std")
 
 # 5 timeslices for PCA
 groups_earlyK <- groupings_help %>% filter(epoch == "Early Cretaceous") 
 groups_midJ <- groupings_help %>% filter(epoch == "Early Cretaceous")
 groups_lateJ <- groupings_help %>% filter(epoch == "Early Cretaceous")
 # groups_JK <- groupings_help %>% filter(stage == "Kimmeridgian,Tithonian, Berriasian, Valanginian")
 # but 13MA vs 9MA
 # last group? midK or earlyJ?


 
 
 
 ####### NEW WAY for PCA
 ## Change row names to occurrence numbers
 PCA_data_earlyK_df <- as.data.frame(PCA_data_earlyK) # switch to dataframe
 rownames(PCA_data_earlyK_df) <- PCA_data_earlyK_df[,1] # copy occurrence_no to row name
 PCA_data_earlyK_df$occurrence_no <- NULL # remove "occurrence_no" column
 head(PCA_data_earlyK_df) #check
 
 
 ## Pull out columns
 earlyK_pca <- PCA_data_earlyK_df[, 1:4]
 earlyK_groups <- as.factor(PCA_data_earlyK_df$family)
 
 
 ## PCA
 earlyK_pca_out <- prcomp(earlyK_pca, center = TRUE, scale. = TRUE) 
 summary(earlyK_pca_out)
 
 earlyK_pca.var <- earlyK_pca_out$sdev^2 # amount of variation each PC counts for
 earlyK_pca.var.per <- round(earlyK_pca.var/sum(earlyK_pca.var)*100, 1) # percentages
 
 barplot(earlyK_pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")
 
 
 ## get the name of the top 10 measurements that contribute most to pc1.
 loading_scores <- earlyK_pca_out$rotation[,1]
 climate_scores <- abs(loading_scores) ## get the magnitudes
 climate_score_ranked <- sort(climate_scores, decreasing=TRUE)
 top_4 <- names(climate_score_ranked[1:4])
 
 earlyK_pca_out$rotation[top_4,1] ## show the scores (and +/- sign)
 
 ## Plot
 earlyK_PCA_plot <- ggbiplot::ggbiplot(earlyK_pca_out, obs.scale = 1, var.scale = 1, groups = earlyK_groups, ellipse = TRUE)
 earlyK_PCA_plot <- earlyK_PCA_plot + scale_colour_manual(values = c("#0891A3", "#FFA93D", "#B00B69", "grey50", "#248528","#572AAC" )) +
   theme(panel.background = element_blank(),
         legend.position = "top", #legend.position = "none",
         panel.border = element_rect(colour = "black", fill = NA))
 earlyK_PCA_plot

