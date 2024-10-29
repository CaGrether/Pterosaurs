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

#### USE ptero_tree_dated from phylo_plots.R !!

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


########## DO I NEED THIS OR HAVE I ALREADY DONE IT WITH THE CODE ABOVE
## Grab a list of species to subset the large dataset
tree_eff_species <- ptero_tree_eff$tip.label

## Remove taxa that are in this dataset but not on the tree:
taxa_to_remove <- species_climate$accepted_name[ !species_climate$accepted_name %in% tree_species ] # in MAT data but not on tree
species_climate_tree <- species_climate[!species_climate$accepted_name %in% taxa_to_remove , ]
species_climate_tree <- species_climate_tree[!species_climate_tree$accepted_name =="Bakonydraco_galaczi", ]

## drop tips for taxa on tree that do not have climate data (there should not be (m)any!)
tree_pruned <- drop.tip(ptero_tree_dated, ptero_tree_dated$tip.label[!(ptero_tree_dated$tip.label %in% species_climate_tree$accepted_name)])
setdiff(tree_pruned$tip.label, species_climate_tree$accepted_name) # taxa on tree, but not in data - should = "character(0)" 

## Get averages of the climate variable for all species:
climate_mean <- species_climate_tree %>% 
  group_by(accepted_name) %>% 
  summarise(mean_MAT = mean(MAT), mean_MAP = mean(MAP))
####################

## turn the species column into the row names
eff_set <- eff_tree %>% 
  group_by(species) 
eff_set <- column_to_rownames(eff_set, var = "species")
eff_set[,10] <- as.numeric(eff_set[,10])


## subset efficiency data for desired variable inv CoT
Inv_CoT_matrix <- as.matrix(eff_set) [,10]
Inv_CoT_matrix <- sapply(Inv_CoT_matrix, as.numeric)

## Efficiency contMap() inv CoT
CoTmapped <- contMap(ptero_tree_eff, Inv_CoT_matrix, plot = FALSE)
CoTmapped <- setMap(CoTmapped, invert = TRUE)
n <- length(CoTmapped$cols)
CoTmapped$cols[1:n] <- plasma(n)
plot(CoTmapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "CoT")




## Precipitation contMap()
MAPmapped <- contMap(tree_pruned, MAP_matrix, plot = FALSE)
MAPmapped <- setMap(MAPmapped, invert = TRUE)
MAPmapped <- drop.tip(MAPmapped, "Bakonydraco_galaczi")
MAPmapped <- drop.tip(MAPmapped, "Eurazhdarcho_langendorfensis")
n <- length(MAPmapped$cols)
MAPmapped$cols[1:n] <- turbo(n, direction = -1)
plot(MAPmapped, fsize = c(0.3, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)", par(bg="#ECDDBF"))


