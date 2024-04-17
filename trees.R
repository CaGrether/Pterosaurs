#### Pterosaurs on PBDB and in Yu et al. 2023


library(phytools)
library(tidyverse)


# read in tree(s) from Yu et al. 2023
cis <- ape::read.nexus("Trees/Data_S1.nex") # Complete informal super tree
cis
psc <- ape::read.nexus("Trees/Data_S2.nex") # parsimony strict consensus tree
psc
pis <- ape::read.nexus("Trees/Data_S3.nex") # pruned informal super tree
pis

# plot
plotTree(cis)
plot(cis, type = "fan", show.tip.label = FALSE)
plot(psc, type = "fan", show.tip.label = FALSE)
plot(pis, type = "fan", show.tip.label = FALSE)

# tip labels (cis from here)
species.yu <- cis$tip.label


## PBDB
# read in PBDB pterosauria
dat <- read.csv("pbdb_data.csv", skip = 18)

# species rank
dat.species <- dplyr::filter(dat, accepted_rank == "species") # only accepted species

# adjust species names to tree
species.pb <- c() # vector for taxa in the tree

for(i in dat.species$accepted_name){
  
  tmp <- gsub(" ", "_", i) # '_' in btw genus and species name for every tip label
  
  # taxa in tree
  species.pb <- c(species.pb, 
                 (ifelse (tmp %in% species.yu, TRUE, FALSE)) )
  
  # cat("\r", i, "of", 100) How to count?
  
}

species.pb

#dat.species$tree_occurrence <- species.pb

# match species.pb and species.yu