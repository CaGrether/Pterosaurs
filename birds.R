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
#   Plots of bird phylogenetic tree from Wang et al. 2021
# 
# ******************************************************

## Load packages
library(ape)
library(phytools)
library(tidyverse)
library(strap)
library(viridis)

# 1. Prep data ------------------------------------------------------------


bird_tree_raw <- read.nexus("Trees/tipdated.tre")

plot(bird_tree_raw, type = "fan", show.tip.label = F)

# get list of taxa
bird_taxa <- bird_tree_raw$tip.label
