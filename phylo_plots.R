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
# 
# ******************************************************


library(ape)
library(phytools)
library(tidyverse)


# ptero_tree <- read.nexus("Trees/Data_S1.nex")
ptero_tree <- read.nexus("Trees/Data_S3.nex") # if using the most parsimonous tree

# get list of taxa
ptero_taxa <- ptero_tree$tip.label
ptero_taxa <- as.data.frame(ptero_taxa)

write_csv(ptero_taxa, "Data/Output/ptero_taxa_s3.csv")



## Get min and max ages for all taxa

## Import the cleaned data (species and ody fossils only)
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20) # your file path will be different!
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

write_csv(ages_info, "Data/Output/ages_info.csv")


## remove specimens
# can't think of a way to do this, so doing it by force
ptero_taxa_clean <- as.data.frame(ptero_taxa[-c(which(ptero_taxa == "OCP_DEK_GE_716"),
                                          which(ptero_taxa == "LPM_L112113"),
                                          which(ptero_taxa == "LPM_N081607")),])
names(ptero_taxa_clean) <- "ptero_taxa"


## get families from PBDB
family_info <- select(occurrences_sp, accepted_name, family)


# merge?
# add underscore
ptero_taxa_vec <- c()

for(i in family_info$accepted_name){ # checking species names in pbdb
  
  tmp <- gsub(" ", "_", i) # change name for every tip label
  
  ptero_taxa_vec <- c(ptero_taxa_vec, tmp)
}

family_info <- as.data.frame(cbind(ptero_taxa_vec, family_info$family))

merge(ptero_taxa_clean, family_info, by.x = "ptero_taxa", 
      by.y = "ptero_taxa_vec", all.x = F)
