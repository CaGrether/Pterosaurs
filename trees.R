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

species.pb # species renamed in PBDB
in.tree # list of T or F whether PBDB are in tree or not

length(unique(species.pb)) # 262 unique species in PBDB
length(which(in.tree == FALSE)) # 74 of 262 not in cis tree
length(unique(species.yu)) # 213 species in tree

# data frame for species in tree
tree.species <- as.data.frame(cbind(species.pb, in.tree))
names(tree.species) <- c("species", "Is in tree")
tree.species

# only FALSE taxa into csv
new.sp <- tree.species[which(tree.species$`Is in tree`==FALSE),]

### I only need the FALSE ones, so only 'new.sp'
new.sp$row_in_pbdb <- rep(NA, length(new.sp$species))
new.sp$ref <- rep(NA, length(new.sp$species))
i <- 1

for (x in new.sp$species) {
  rownr <- as.character(which(dat.species$tree_name == x))
  
  # get row number in original pbdb data
  new.sp$row_in_pbdb[i] <- paste(rownr, sep = "", collapse = ",")
  
    # get reference (author, year) -> ref_author, ref_pubyr
    ath.yr <- c() # author and year together
    
    for (v in length(rownr)) {
      tmp <- as.integer(rownr[v])
    
      basket <- paste(dat.species$ref_author[tmp], dat.species$ref_pubyr[tmp], 
                           sep = "", collapse = ",") # author and year for each specimen
      ath.yr <- as.vector(c(ath.yr, basket)) # author(s) and year(s) for species
    }
  
  # add reference into new species file  
  new.sp$ref[i] <- paste(unique(ath.yr), sep = "", collapse = ",")
  
  # set counter
  i <- i+1
}

# remove FALSE column
new.sp <- new.sp[,-2]

#write.csv2(new.sp, "Data/output/cis_missing.csv")
#write_csv(new.sp, "Data/output/cis_missing2.csv") # for compatibility with Emma's computer lol! 

################################################################################
#### why are they not in tree? Compare with Yu et al supplementary S5 listing
#### missing species and the reason for why they're missing
################################################################################

# compare with online data S5 
s5 <- read.csv2("Data/Input/Data_S5_copy.csv", skip = 1) # s5 data
sp_missing<- read.csv2("Data/Output/cis_missing.csv") # missing species from tree

# unique species in s5
s5_sp <- unique(s5$Name)

length(s5_sp)
length(sp_missing$species)

miss_inS5 <-  # missing species that are (explained) in S5
  sp_missing$species[which(sp_missing$species %in% s5_sp)]  


# only need one reason for why it's not in (?)
# don't need reason for every single specimen

reason_miss <- s5[match(miss_inS5, s5$Name),]
# many with NA, but why Peteinosaurus_zambelli in analysis but not in tree?

#Name <- c()
#for (x in miss_inS5) {
#  Name <- c(Name, rep(x, length(which(s5$Name == x))))
#}


