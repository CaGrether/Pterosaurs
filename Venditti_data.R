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

## Load the tree and data files:
ptero_tree_dated <- read.nexus("Trees/ptero_tree_dated.nex")
eff.dat <- read.csv2("Data/Output/Species_eff_data.csv")

## Drop two outlier taxa Bakonydraco galaczi Eurazhdarcho langendorfensis (if wanted)

# ptero_tree_dated <- drop.tip(ptero_tree_dated, "Eurazhdarcho_langendorfensis")
# ptero_tree_dated <- drop.tip(ptero_tree_dated, "Bakonydraco_galaczi")

#### analysis
# To include Quetzalcoatlus data, change eff.dat "Quetzalcoatlus spp" to "Quetzalcoatlus northropi"
eff.dat$species[62] <- "Quetzalcoatlus_northropi"

eff.dat$species[36] <- "Hatzegopteryx_thambema" # fixed spelling
eff.dat$species[40] <- "Huaxiapterus_corollatus" # fixed spelling

# compare species in tree and Venditti data
setdiff(eff.dat$species, ptero_tree_dated$tip.label) # taxa in data but not on tree 
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
#CoTmapped <- drop.tip(CoTmapped, "Bakonydraco_galaczi")
#CoTmapped <- drop.tip(CoTmapped, "Eurazhdarcho_langendorfensis")
n <- length(CoTmapped$cols)
CoTmapped$cols[1:n] <- rocket(n)
plot(CoTmapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "CoT (kg m J")


# ______________________________________________________
#
#   Phylogenetic signal in flight efficiency data (Venditti et al. 2020)
# 
# ******************************************************

## Test the phlogenetic signal (Pagel's λ and Blomberg's K)
lambdaEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "lambda", test = TRUE)
lambdaEff 
# P-Value >0.001 and Lambda 0.9999

KEff <- phylosig(ptero_tree_eff, Inv_CoT_matrix, method = "K", test = TRUE, nsim = 10000)
KEff
# P-Value >0.001 and K 2.58778

##### exclude pteranodontoidea group
tip <- c("Pteranodon_sternbergi", "Anhanguera_blittersdorffi")
test.tree <- drop.clade(ptero_tree_eff, tip)
test.tree <- drop.tip(test.tree, "NA")

## remove species from efficiency data
spec_to_remove <- setdiff(ptero_tree_eff$tip.label, test.tree$tip.label) 

tmp <- rownames_to_column(eff_set, var = "species")
tmp <- tmp[!tmp$species %in% spec_to_remove, ]

eff_set_pruned <- tmp %>% 
          group_by(species) %>% 
          column_to_rownames(var = "species")

setdiff(test.tree$tip.label, tmp$species)

## subset efficiency data for desired variable inv CoT
Inv_CoT_pruned <- as.matrix(eff_set_pruned) [,10]
Inv_CoT_pruned <- sapply(Inv_CoT_pruned, as.numeric)

## Test the phlogenetic signal again
lambdaEffpruned <- phylosig(test.tree, Inv_CoT_pruned, method = "lambda", test = TRUE)
lambdaEffpruned

KEffpruned <- phylosig(test.tree, Inv_CoT_pruned, method = "K", test = TRUE, nsim = 10000)
KEffpruned

# plot tree without Pteranodontoidea
## Efficiency contMap() inv CoT
CoTpruned <- contMap(test.tree, Inv_CoT_pruned, plot = FALSE)
CoTpruned <- setMap(CoTpruned, invert = TRUE)
#CoTpruned <- drop.tip(CoTpruned, "Bakonydraco_galaczi")
#CoTpruned <- drop.tip(CoTpruned, "Eurazhdarcho_langendorfensis")
n <- length(CoTpruned$cols)
CoTpruned$cols[1:n] <- plasma(n)
plot(CoTpruned, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "CoT")


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
#species_climEff_tree <- species_climEff_tree[!species_climEff_tree$accepted_name =="Bakonydraco_galaczi", ]
#species_climEff_tree <- species_climEff_tree[!species_climEff_tree$accepted_name =="Eurazhdarcho_langendorfensis", ]

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
MATmapped_eff$cols[1:n] <- plasma(n)
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
Tseason_mapped$cols[1:n] <- plasma(n)
plot(Tseason_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "T (°C)")

## PSeason contMap()
Pseason_mapped <- contMap(tree_eff_pruned, Pseason_matrix, plot = FALSE)
Pseason_mapped <- setMap(Pseason_mapped, invert = TRUE)
n <- length(Pseason_mapped$cols)
Pseason_mapped$cols[1:n] <- viridis(n, direction = -1)
plot(Pseason_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "P (mm/day)")


### remove extreme Cuban and German species
# get these species
xtaxa <- c("Nesodactylus_hesperius", "Rhamphorhynchus_muensteri")
xf.tree <- drop.clade(tree_eff_pruned, xtaxa)
xf.tree <- drop.tip(xf.tree, "NA")

# get list of taxa in that tree
xftree_spec <- xf.tree$tip.label

## Remove taxa that are in climate data but not on the tree:
taxa_to_remove <- species_climate$accepted_name[ !species_climate$accepted_name %in% xftree_spec ] # in MAT data but not on tree
species_climxf_tree <- species_climate[!species_climate$accepted_name %in% taxa_to_remove , ]

## drop tips for taxa on tree that do not have climate data (there should not be (m)any!)
tree_xf_pruned <- drop.tip(xf.tree, xf.tree$tip.label[!(xf.tree$tip.label %in% species_climxf_tree$accepted_name)])
setdiff(tree_xf_pruned$tip.label, species_climxf_tree$accepted_name) # taxa on tree, but not in data - should = "character(0)" 


## Get averages of the climate variable for all species:
climate_xf_mean <- species_climxf_tree %>% 
  group_by(accepted_name) %>% 
  summarise(mean_MAT = mean(MAT), mean_MAP = mean(MAP), 
            mean_Tseason = mean(seasonal_temp), mean_Pseason = mean(seasonal_precip))

## turn the accepted_name column into the row names
climate_xf_mean <- column_to_rownames(climate_xf_mean, var = "accepted_name")

## Convert to matrix:
MAT_xf_matrix <- as.matrix(climate_xf_mean) [,1] # Mean annual temperature
MAP_xf_matrix <- as.matrix(climate_xf_mean) [,2] # Mean annual precipitation
Tseason_xf_matrix <- as.matrix(climate_xf_mean) [,3] # Mean seasonal temperature
Pseason_xf_matrix <- as.matrix(climate_xf_mean) [,4] # Mean seasonal precipitation

## Temperature contMap()
MATmapped_xf <- contMap(tree_xf_pruned, MAT_xf_matrix, plot = FALSE)
MATmapped_xf <- setMap(MATmapped_xf, invert = TRUE)
n <- length(MATmapped_xf$cols)
MATmapped_xf$cols[1:n] <- plasma(n)
plot(MATmapped_xf, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAT (°C)")

## Precipitation contMap()
MAPmapped_xf <- contMap(tree_xf_pruned, MAP_xf_matrix, plot = FALSE)
MAPmapped_xf <- setMap(MAPmapped_xf, invert = TRUE)
n <- length(MAPmapped_xf$cols)
MAPmapped_xf$cols[1:n] <- viridis(n, direction = -1)
plot(MAPmapped_xf, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "MAP (mm/day)")

## TSeason contMap()
Tseason_xf_mapped <- contMap(tree_xf_pruned, Tseason_xf_matrix, plot = FALSE)
Tseason_xf_mapped <- setMap(Tseason_xf_mapped, invert = TRUE)
n <- length(Tseason_xf_mapped$cols)
Tseason_xf_mapped$cols[1:n] <- plasma(n)
plot(Tseason_xf_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "T (°C)")

## PSeason contMap()
Pseason_xf_mapped <- contMap(tree_xf_pruned, Pseason_xf_matrix, plot = FALSE)
Pseason_xf_mapped <- setMap(Pseason_xf_mapped, invert = TRUE)
n <- length(Pseason_xf_mapped$cols)
Pseason_xf_mapped$cols[1:n] <- viridis(n, direction = -1)
plot(Pseason_xf_mapped, fsize = c(0.4, 1), outline = FALSE, lwd = c(3, 7), leg.txt = "P (mm/day)")
