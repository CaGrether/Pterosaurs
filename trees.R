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
# 1. for new names in data frame
#### need column with the accepted name written with underscore, not unique but for all 440
#### so that I can get all references into one column in my csv Not_in_cis
#### problem: need to do 2nd for loop to make it not unique
#### but now need to go through names but also know the rownumber for each name how??
treename <- c()

for(i in dat.species$accepted_name){ # checking species names in pbdb
  
  tmp <- gsub(" ", "_", i) # change name for every tip label
  
  treename <- c(treename, tmp)
}

dat.species$tree_name <- treename # put into pbdb data for later reference

# 2. for actual analysis
species.pb <- c() # vector for taxa in the tree
in.tree <- c() # vector for T/F

for(i in unique(dat.species$accepted_name)){ # checking unique species names in pbdb
  
  tmp <- gsub(" ", "_", i) # change name for every tip label
  
  species.pb <- c(species.pb, tmp) # store new name
  
  # is species in Yu et al. tree
  in.tree <- c(in.tree, 
                 (ifelse (tmp %in% species.yu, TRUE, FALSE)) )
  
  # cat("\r", i, "of", 100) How to count?
  
}

species.pb
in.tree

length(unique(species.pb)) # 262 unique species in PBDB
length(which(in.tree == FALSE)) # 74 of 262 not in cis tree
length(unique(species.yu)) # 213 species in tree

# data frame for species in tree
tree.species <- as.data.frame(cbind(species.pb, in.tree))
names(tree.species) <- c("species", "Is in tree")
tree.species

# only FALSE taxa into csv
new.sp <- tree.species[which(tree.species$`Is in tree`==FALSE),]
#write.csv2(new.sp, "Data/Not_in_cis_tree.csv")


### I only need the FALSE ones, so only 'new.sp'
test <- c()

for (x in new.sp$species) {
  rownr <- which(dat.species$tree_name == x)
  test <- c(test, rownr) # somehow differentiate between multiple occurrences
}

# later
ref <- data.frame(species.pb, rep(NA, length(species.pb)))