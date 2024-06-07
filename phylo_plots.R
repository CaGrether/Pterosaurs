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


ptero_tree <- read.nexus("Trees/Data_S1.nex")

# get list of taxa
ptero_taxa <- ptero_tree$tip.label
ptero_taxa <- as.data.frame(ptero_taxa)

write_csv(ptero_taxa, "Data/Output/ptero_taxa.csv")



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


## get stratoographic ranges:
ages_info <- taxon_info %>% group_by(accepted_name) %>% 
  summarise(FAD = max(max_ma), LAD = min(min_ma))





