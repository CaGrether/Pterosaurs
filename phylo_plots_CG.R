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
#   Plots of climate on phylogenetic tree
#   full species tree
# 
# ******************************************************

## Load packages
library(ape)
library(phytools)
library(tidyverse)
library(strap)
library(viridis)
library(desk) # for Jarque Bera test



# 1. Prep data ------------------------------------------------------------


ptero_tree_raw <- read.nexus("Trees/Data_S3.nex")

# get list of taxa
ptero_taxa <- ptero_tree_raw$tip.label
ptero_taxa <- as.data.frame(ptero_taxa)

write_csv(ptero_taxa, "Data/Output/ptero_taxa_s3.csv")


## Get min and max ages for all taxa

## Import the cleaned data (species and ody fossils only)
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20) # fix file path if necessary!
glimpse(occurrences)

## species only:
occurrences_sp <- occurrences %>% filter(accepted_rank == "species")


## take necessary columns
taxon_info <- select(occurrences_sp, collection_no, occurrence_no, accepted_name, 
                       early_interval, late_interval, max_ma, min_ma)
glimpse(taxon_info)


## get stratographic ranges:
ages_info <- taxon_info %>% group_by(accepted_name) %>% 
  summarise(FAD = max(max_ma), LAD = min(min_ma))
glimpse(ages_info)

## save a copy:
write_csv(ages_info, "Data/Output/ages_info.csv")



# 2. Load & organise tree files ------------------------------------------------

## Load the tree file:
ptero_tree_raw <- read.nexus("Trees/Data_S3.nex")

## Replace spaces with underscores in the age data (see above) to match the tree
ages_info$accepted_name <- gsub(" ", "_", ages_info$accepted_name)

## Check that the tree and the age data match
setdiff(ages_info$accepted_name, ptero_tree_raw$tip.label) # taxa in data but not on tree
setdiff(ptero_tree_raw$tip.label, ages_info$accepted_name) # taxa on tree, but not in data

## Drop tips from tree that have no data in the dataset
ptero_tree <- drop.tip(ptero_tree_raw, ptero_tree_raw$tip.label[!(ptero_tree_raw$tip.label %in% ages_info$accepted_name)])
setdiff(ptero_tree$tip.label, ages_info$accepted_name) # taxa on tree, but not in data - should = "character(0)" 

## Drop tips from data that are not on the tree
ages_tree <- ages_info[(ages_info$accepted_name %in% ptero_tree_raw$tip.label), ]
setdiff(ages_tree$accepted_name, ptero_tree_raw$tip.label) # taxa in data but not on tree - should = "character(0)" 

## move accepted_name column to rownames
ptero_timeData <- ages_tree %>% column_to_rownames(var = "accepted_name")

## Add the ages to the tree to 'date' it
ptero_tree_dated <- DatePhylo(ptero_tree, ptero_timeData, method = "equal", rlen = 1)

## Drop two outlier taxa Bakonydraco galaczi Eurazhdarcho langendorfensis

#ptero_tree_dated <- drop.tip(ptero_tree_dated, "Eurazhdarcho_langendorfensis")
#ptero_tree_dated <- drop.tip(ptero_tree_dated, "Bakonydraco_galaczi")

## save tree for later purpose
writeNexus(ptero_tree_dated, file = "Trees/ptero_tree_dated.nex")

# 3. Organise the climate data --------------------------------------------

## Import climate data:
species_climate <- read_csv("Data/climate/species_climate.csv")

## Replace spaces with underscores in the species climate data (see above) to match the tree
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)


## Grab a list of species to subset the large dataset
tree_species <- ptero_tree$tip.label

## Remove taxa that are in this dataset but not on the tree:
taxa_to_remove <- species_climate$accepted_name[ !species_climate$accepted_name %in% tree_species ] # in MAT data but not on tree
species_climate_tree <- species_climate[!species_climate$accepted_name %in% taxa_to_remove , ]

## drop tips for taxa on tree that do not have climate data (there should not be (m)any!)
tree_pruned <- drop.tip(ptero_tree_dated, ptero_tree_dated$tip.label[!(ptero_tree_dated$tip.label %in% species_climate_tree$accepted_name)])
setdiff(tree_pruned$tip.label, species_climate_tree$accepted_name) # taxa on tree, but not in data - should = "character(0)" 


## Get averages of the climate variable for all species:
climate_mean <- species_climate_tree %>% 
  group_by(accepted_name) %>% 
  summarise(mean_MAT = mean(MAT), mean_MAP = mean(MAP), 
            mean_T = mean(seasonal_temp), mean_P = mean(seasonal_precip))

## turn the accepted_name column into the row names
climate_mean <- column_to_rownames(climate_mean, var = "accepted_name")

## Convert to matrix:
MAT_matrix <- as.matrix(climate_mean) [,1] # Mean annual temperature
MAP_matrix <- as.matrix(climate_mean) [,2] # Mean annual precipitation
T_matrix <-as.matrix(climate_mean) [,3] # Mean seasonal temperature
P_matrix <- as.matrix(climate_mean) [,4] # Mean seasonal precipitation


# Jarque Bera test for normal distribution
jbMAT <- jb.test(MAT_matrix)
jbMAP <- jb.test(MAP_matrix)
jbT <- jb.test(T_matrix)
jbP <- jb.test(P_matrix)

## Temperature contMap()
MATmapped <- contMap(tree_pruned, MAT_matrix, plot = FALSE)
MATmapped <- setMap(MATmapped, invert = TRUE)
#MATmapped <- drop.tip(MATmapped, "Bakonydraco_galaczi")
#MATmapped <- drop.tip(MATmapped, "Eurazhdarcho_langendorfensis")
n <- length(MATmapped$cols)
MATmapped$cols[1:n] <- plasma(n)
plot(MATmapped, fsize = c(0.4, 1), fcol = "red", outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")


## Precipitation contMap()
MAPmapped <- contMap(tree_pruned, MAP_matrix, plot = FALSE)
MAPmapped <- setMap(MAPmapped, invert = TRUE)
#MAPmapped <- drop.tip(MAPmapped, "Bakonydraco_galaczi")
#MAPmapped <- drop.tip(MAPmapped, "Eurazhdarcho_langendorfensis")
n <- length(MAPmapped$cols)
MAPmapped$cols[1:n] <- viridis(n, direction = -1)
plot(MAPmapped, fsize = c(0.3, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)")


## seasonal Temperature contMap
Tmapped <- contMap(tree_pruned, T_matrix, plot = FALSE)
Tmapped <- setMap(Tmapped, invert = TRUE)
#Tmapped <- drop.tip(MAPmapped, "Bakonydraco_galaczi")
#Tmapped <- drop.tip(MAPmapped, "Eurazhdarcho_langendorfensis")
n <- length(Tmapped$cols)
Tmapped$cols[1:n] <- plasma(n)
plot(Tmapped, fsize = c(0.5, 1), fcol = "red", outline = FALSE, lwd = c(3, 7), leg.txt = "T (°C)")


## seasonal Precipitation contMap
Pmapped <- contMap(tree_pruned, P_matrix, plot = FALSE)
Pmapped <- setMap(Pmapped, invert = TRUE)
#Pmapped <- drop.tip(MAPmapped, "Bakonydraco_galaczi")
#Pmapped <- drop.tip(MAPmapped, "Eurazhdarcho_langendorfensis")
n <- length(Pmapped$cols)
Pmapped$cols[1:n] <- viridis(n, direction = -1)
plot(Pmapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "P (mm/day)")


# phylogenetic signal
lambdaMAT <- phylosig(tree_pruned, MAT_matrix, method = "lambda", test = TRUE)
lambdaMAT 

lambdaMAP <- phylosig(tree_pruned, MAP_matrix, method = "lambda", test = TRUE)
lambdaMAP 

lambdaT <- phylosig(tree_pruned, T_matrix, method = "lambda", test = TRUE)
lambdaT 

lambdaP <- phylosig(tree_pruned, P_matrix, method = "lambda", test = TRUE)
lambdaP 


KMAT <- phylosig(tree_pruned, MAT_matrix, method = "K", test = TRUE, nsim = 10000)
KMAT 

KMAP <- phylosig(tree_pruned, MAP_matrix, method = "K", test = TRUE, nsim = 10000)
KMAP 

KT <- phylosig(tree_pruned, T_matrix, method = "K", test = TRUE, nsim = 10000)
KT 

KP <- phylosig(tree_pruned, P_matrix, method = "K", test = TRUE, nsim = 10000)
KP 


##### exclude extreme group of Cuban/German pterosaurs
xtaxa <- c("Nesodactylus_hesperius", "Rhamphorhynchus_muensteri")
x.tree <- drop.clade(tree_pruned, xtaxa)
x.tree <- drop.tip(x.tree, "NA")

# get list of taxa in that tree
xtree_spec <- x.tree$tip.label

# fit climate data to new number of taxa
## Remove taxa that are in climate data but not on the tree:
taxa_to_remove <- species_climate$accepted_name[ !species_climate$accepted_name %in% xtree_spec ] # in data but not on tree
xspec_climatree <- species_climate[!species_climate$accepted_name %in% taxa_to_remove , ]

## Get averages of the climate variable for xtree species:
xclimate_mean <- xspec_climatree %>% 
  group_by(accepted_name) %>% 
  summarise(mean_MAT = mean(MAT), mean_MAP = mean(MAP), 
            mean_T = mean(seasonal_temp), mean_P = mean(seasonal_precip))

## turn the accepted_name column into the row names
xclimate_mean <- column_to_rownames(xclimate_mean, var = "accepted_name")

## Convert to matrix:
MATx_matrix <- as.matrix(xclimate_mean) [,1]
MAPx_matrix <- as.matrix(xclimate_mean) [,2]
Tx_matrix <- as.matrix(xclimate_mean) [,3]
Px_matrix <- as.matrix(xclimate_mean) [,4] # seasonal Preciptiation

## MAT (see impact)
MAT_x_mapped <- contMap(x.tree, MATx_matrix, plot = FALSE)
MAT_x_mapped <- setMap(MAT_x_mapped, invert = TRUE)
n <- length(MAT_x_mapped$cols)
MAT_x_mapped$cols[1:n] <- plasma(n)
plot(MAT_x_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")

## MAP (see impact)
MAP_x_mapped <- contMap(x.tree, MAPx_matrix, plot = FALSE)
MAP_x_mapped <- setMap(MAP_x_mapped, invert = TRUE)
n <- length(MAP_x_mapped$cols)
MAP_x_mapped$cols[1:n] <- viridis(n, direction = -1)
plot(MAP_x_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)")

## seasonal Temperature (see impact)
T_x_mapped <- contMap(x.tree, Tx_matrix, plot = FALSE)
T_x_mapped <- setMap(T_x_mapped, invert = TRUE)
n <- length(T_x_mapped$cols)
T_x_mapped$cols[1:n] <- plasma(n)
plot(T_x_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "T (°C)")

## seasonal Precipitation (see impact)
P_x_mapped <- contMap(x.tree, Px_matrix, plot = FALSE)
P_x_mapped <- setMap(P_x_mapped, invert = TRUE)
n <- length(P_x_mapped$cols)
P_x_mapped$cols[1:n] <- viridis(n, direction = -1)
plot(P_x_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "P (mm/day)")




# 4. Select groups of pterosaurs (work in progress) ----------------------------------------

# ## remove specimens
# ptero_taxa_clean <- as.data.frame(ptero_taxa[-c(which(ptero_taxa == "OCP_DEK_GE_716"),
#                                           which(ptero_taxa == "LPM_L112113"),
#                                           which(ptero_taxa == "LPM_N081607")),])
# names(ptero_taxa_clean) <- "ptero_taxa"
# 
# 
# ## get families from PBDB
# family_info <- select(occurrences_sp, accepted_name, family)
# 
# 
# ## add underscore to match names
# family_info$accepted_name <- gsub(" ", "_", family_info$accepted_name)
# 
# 
# ## merge
# ptero_fam <- merge(ptero_taxa_clean, family_info, by.x = "ptero_taxa", 
#       by.y = "accepted_name")
# 
# 
# ## remove duplicates
# ptero_groups <- ptero_fam[!duplicated(ptero_fam),]
# names(ptero_groups) <- c("ptero_taxa", "family")
# 
# 
# ## save a copy
# write_csv(ptero_groups, "Data/Output/ptero_groups.csv")
