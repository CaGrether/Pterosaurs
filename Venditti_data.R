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
#   Venditti et al. 2020 paper data prep
#   Contmapping
# 
# ******************************************************

########## data preparation ----------------------------------------------

# load in data
eff.dat <- read.table("Data/Input/Venditti_avgData.txt")

# new header
names(eff.dat)
names(eff.dat) <- eff.dat[1,]
names(eff.dat)

eff.dat <- eff.dat[-1,]

write.csv2(eff.dat,"Data/Output/Species_eff_data.csv", row.names = FALSE)


########## analysis ------------------------------------------------------

library(phytools)
library(tidyverse)
library(strap)
library(viridis)
library(ape)


# 1. Load & organise tree file ------------------------------------------------

## Load the tree file:
ptero_tree_raw <- read.nexus("Trees/Data_S3.nex")
ages_info <- read.csv("Data/Output/ages_info.csv")

#### Use ptero_tree_dated from phylo_plots.R 

#### phylo_plots.R Script modified
# compare species in tree and Venditti data
setdiff(eff.dat$species, ptero_tree_dated$tip.label) # taxa in data but not on tree 1 (Quetzalcoatlus spp)
setdiff(ptero_tree_dated$tip.label, eff.dat$species) # taxa on tree, but not in data 82

# Drop tips from data that are not on the tree
eff_tree <- eff.dat[(eff.dat$species %in% ptero_tree_dated$tip.label), ]
setdiff(eff_tree$species, ptero_tree_dated$tip.label) # taxa in data but not on tree - should = "character(0)" 

# Drop tips from tree that are not in the data
ptero_tree_eff <- drop.tip(ptero_tree_dated, ptero_tree_dated$tip.label[!(ptero_tree_dated$tip.label %in% eff.dat$species)])
setdiff(eff_tree$species, ptero_tree_dated$tip.label) # taxa in data but not on tree - should = "character(0)" 


## turn the species column into the row names
eff_set <- eff_tree %>% 
  group_by(species) 
eff_set <- column_to_rownames(eff_set, var = "species")

## subset efficiency data for desired variable inv CoT
Inv_CoT_matrix <- as.matrix(eff_set) [,10]
Inv_CoT_matrix <- sapply(Inv_CoT_matrix, as.numeric)

## Efficiency contMap() inv CoT
CoTmapped <- contMap(ptero_tree_eff, Inv_CoT_matrix, plot = FALSE)
CoTmapped <- setMap(CoTmapped, invert = TRUE)
n <- length(CoTmapped$cols)
CoTmapped$cols[1:n] <- plasma(n)
plot(CoTmapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "CoT")


# ______________________________________________________
#
#   Phylogenetic signal in flight efficiency data (Venditti et al. 2020)
# 
# ******************************************************

# Prepare for analysis

library(phytools)
library(tidyverse)
library(strap)
library(viridis)
library(ape)

## Test the phlogenetic signal (Pagel's λ and Blomberg's K)
lambdaEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "lambda", test = TRUE)
lambdaEff 
# P-Value >0.001 and Lambda 0.9999

KEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "K", test = TRUE, nsim = 10000)
KEff
# P-Value >0.001 and K 2.58778


# Combine with climate data
# ______________________________________________________
#
#   Plots of climate on Venditti tree
#   efficiency data - pruned tree
# 
# ******************************************************

## Import climate data:
species_climate <- read_csv("Data/climate/species_climate.csv")

## Replace spaces with underscores in the species climate data (see above) to match the tree
species_climate$accepted_name <- gsub(" ", "_", species_climate$accepted_name)

## Grab a list of species to subset the large dataset
tree_eff_species <- ptero_tree_eff$tip.label

## Remove taxa that are in this dataset but not on the tree:
taxa_to_remove <- species_climate$accepted_name[ !species_climate$accepted_name %in% tree_eff_species ] # in MAT data but not on tree
species_climEff_tree <- species_climate[!species_climate$accepted_name %in% taxa_to_remove , ]
species_climEff_tree <- species_climEff_tree[!species_climEff_tree$accepted_name =="Bakonydraco_galaczi", ]
species_climEff_tree <- species_climEff_tree[!species_climEff_tree$accepted_name =="Eurazhdarcho_langendorfensis", ]

## drop tips for taxa on tree that do not have climate data (there should not be (m)any!)
tree_eff_pruned <- drop.tip(ptero_tree_eff, ptero_tree_eff$tip.label[!(ptero_tree_eff$tip.label %in% species_climEff_tree$accepted_name)])
setdiff(tree_eff_pruned$tip.label, species_climEff_tree$accepted_name) # taxa on tree, but not in data - should = "character(0)" 


## Get averages of the climate variable for all species:
climate_eff_mean <- species_climEff_tree %>% 
  group_by(accepted_name) %>% 
  summarise(mean_MAT = mean(MAT), mean_MAP = mean(MAP), 
            mean_Tseason = mean(seasonal_temp), mean_Pseason = mean(seasonal_precip))

## turn the accepted_name column into the row names
climate_eff_mean <- column_to_rownames(climate_eff_mean, var = "accepted_name")

## Convert to matrix:
MAT_eff_matrix <- as.matrix(climate_eff_mean) [,1] # Mean annual temperature
MAP_eff_matrix <- as.matrix(climate_eff_mean) [,2] # Mean annual precipitation
Tseason_matrix <- as.matrix(climate_eff_mean) [,3] # Mean seasonal temperature
Pseason_matrix <- as.matrix(climate_eff_mean) [,4] # Mean seasonal precipitation

## Temperature contMap()
MATmapped_eff <- contMap(tree_eff_pruned, MAT_eff_matrix, plot = FALSE)
MATmapped_eff <- setMap(MATmapped_eff, invert = TRUE)
n <- length(MATmapped_eff$cols)
MATmapped_eff$cols[1:n] <- inferno(n)
plot(MATmapped_eff, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")

## Precipitation contMap()
MAPmapped_eff <- contMap(tree_eff_pruned, MAP_eff_matrix, plot = FALSE)
MAPmapped_eff <- setMap(MAPmapped_eff, invert = TRUE)
n <- length(MAPmapped_eff$cols)
MAPmapped_eff$cols[1:n] <- viridis(n, direction = -1)
plot(MAPmapped_eff, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)")

## TSeason contMap()
Tseason_mapped <- contMap(tree_eff_pruned, Tseason_matrix, plot = FALSE)
Tseason_mapped <- setMap(Tseason_mapped, invert = TRUE)
n <- length(Tseason_mapped$cols)
Tseason_mapped$cols[1:n] <- inferno(n)
plot(Tseason_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")

## PSeason contMap()
Pseason_mapped <- contMap(tree_eff_pruned, Pseason_matrix, plot = FALSE)
Pseason_mapped <- setMap(Pseason_mapped, invert = TRUE)
n <- length(Pseason_mapped$cols)
Pseason_mapped$cols[1:n] <- viridis(n, direction = -1)
plot(Pseason_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)")
