# ******************************************************
#
#   Master Thesis
#
#   Occurrences of pterosaurs
#
#   Carolin Grether
# ______________________________________________________
#
#   Various statistics on raw occurrence data
# 
# ******************************************************


# 0. Packages used in this script -----------------------------------------

library(tidyverse)
library(geoscale) # for plotting with the geological time scale on the x-axis (uses base R syntax)
library(viridis) # for colour scales
library(vegan) # for diversity metrics
library(ggplot2) # for plotting
library(ggpubr) # for plotting
library(divDyn) # stages info
library(dplyr) # functions

data(stages) # age

select <- dplyr::select # ensure the select function is coming from dplyr

# ------------------------------------------------------------
# Data
occurrences <- read_csv("Data/Input/pbdb_pterosauromorpha.csv", skip = 20)
occurrences_sp <- occurrences %>% filter(accepted_rank == "species")
## take necessary columns
taxon_inf <- select(occurrences_sp, collection_no, occurrence_no, accepted_name, 
                    early_interval, late_interval, max_ma, min_ma, formation, 
                    collection_name, lagerstatten)



ints_standard <- read.csv2("Data/Input/ints_standard_copy.csv") # manually fixed data
late_int <- read.csv2("Data/Input/stages_for_climate_data_late.csv")

# 0. put data into correct format -----------------------------------------------------------

## add stage names to data
Adat <- left_join(taxon_inf, ints_standard, by = "early_interval")
Bdat <- left_join(taxon_inf, late_int, by = "late_interval")

# merge datasets
big.dat <- rbind(Adat,Bdat)
# remove duplicates
taxon.dat <- big.dat[!duplicated(big.dat),]

## add ages to intervals
meso <- stages[52:81,]
intervals <- meso[,c(4,6,7,8)]

#  Sampling proxy counts

## calculate counts of sampling proxies and plot these alongside raw diversity

# Taxa per interval 
count_taxa <- vector("numeric") # create empty vector for the loop below to populate
for (i in 1:nrow(intervals)) { # for-loop to count each taxon that appears in each interval
  out <- taxon.dat[which(taxon.dat$stage==intervals$stage[i]),]
  count_taxa[i] <- (length(unique(out$accepted_name)))
  print(count_taxa[i])
}

# Collections per interval
count_colls <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- taxon.dat[which(taxon.dat$stage==intervals$stage[i]),]
  count_colls[i] <- (length(unique(out$collection_no)))
  print(count_colls[i])
}

# Formations per interval
count_formations <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- taxon.dat[which(taxon.dat$stage==intervals$stage[i]),]
  count_formations[i] <- (length(unique(out$formation)))
  print(count_formations[i])
}


## Gather the proxy information together in a new dataframe for plotting:
proxy_counts <- data.frame(intervals$stage, intervals$mid, count_taxa, count_colls, count_formations) 
## Rename the columns
proxy_counts <- rename(proxy_counts, 
                       "interval_name" = "intervals.stage", 
                       "mid_ma" = "intervals.mid")

## Convert all zero's to NAs for plotting 
proxy_counts[proxy_counts == 0] <- NA 



# 1. Sampling plots ----------------------------------------------------------

## Option 1: Plotting using ggplot

## Set interval boundaries for the dotted lines on the plot
int_boundaries <- c( 250, 241, 209, 183, 165, 149, 133, 113, 89, 72)
                    

## Set up ggplot layers
proxy_plot <- ggplot() + 
  # Formations (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", linewidth = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", size = 4, shape = 16) +
  # Collections (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", linewidth = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", size = 5, shape = 16) +
  # Taxa (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_taxa), colour = 'black', linewidth = 1.2)  +
  geom_point(data = proxy_counts, aes(mid_ma, count_taxa), colour = "black", size = 4, shape = 16) +
  # Add a minimal theme
  theme_minimal() + 
  labs(x = "Time (Ma)", y = "Sampling proxy counts") +
  # reverse the x-axis to match geological time
  scale_x_reverse(breaks = int_boundaries) +
  # y-axis with even breaks that match the totals in the dataframe:
  scale_y_continuous(breaks = seq(0, 100, 20))
## Call the finished plot
proxy_plot

## Set dimensions and save plot
ggsave(plot = proxy_plot,
       width = 20, height = 15, dpi = 500, units = "cm", 
       filename = "./plots/sampling_proxies.pdf", useDingbats=FALSE)



## Option 2: Plotting using geoscale

## set up parameters for exporting a PDF of the plot
pdf("./plots/sampling_proxies_geoscale.pdf", width = 9, height = 8) 

## Set up the base of the plot with the timescale:
geoscalePlot(proxy_counts$mid_ma, proxy_counts$count_taxa, # ages and main data points
             units = c("Period", "Age"), # which intervals to show in the timeline  
             tick.scale = "Epoch", # resolution of the tick marks on timescale
             boxes = "Age", # option to include grey boxes for individual time bins
             abbrev = c("Age"), # option to abbreviate names of geological units
             lty = 1, pch = NA, 
             cex.age = 0.6, cex.ts = 1, cex.pt = 1, # size of numbers, text, and points on the scale bar
             age.lim = c(245.5, 72.0), # oldest and youngest ages of the entire time interval
             data.lim = c(0, 100), # range of the data
             ts.col = TRUE, # include colours in the timescale 
             ts.width = 0.17, # space taken up by plotting the time scale
             label = "Counts", # label for y-axis
             direction ="horizontal", # orientation of the plot
             erotate = 0) # numerical value for the rotation for the temporal units

# Add the different lines and points separately:
# Formations (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_formations, col="peru", lty = 2, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_formations, pch =16, cex = 1.5, col = "peru")
# Collections (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_colls, col="orange", lty = 3, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_colls, pch =16, cex = 1.5, col = "orange")
# Taxa (as dots and a line):
lines(proxy_counts$mid_ma, proxy_counts$count_taxa, col="black", lty = 1, lwd = 3)
points(proxy_counts$mid_ma, proxy_counts$count_taxa, pch =16, cex = 1.5, col = "black")

# Finally, add a legend
legend('topright', legend = c("Formations", "Collections", "Taxa"), 
       col = c("peru", "orange", "black"), 
       pch = c(16), box.lty = 1, pt.cex=2, bg= "white")

dev.off() ## Turn off graphic device, to trigger the export of the pdf

### find and compare number of Lagerstätten in Aptian and Campanian 
# filter to stage
tax.dat_Apt <- taxon.dat %>% filter(stage == "Aptian")
tax.dat_Cam <- taxon.dat %>% filter(stage == "Campanian")
tax.dat_Toa <- taxon.dat %>% filter(stage == "Toarcian")

# unique names of Lagerstätten
L_names_Apt <- unique(tax.dat_Apt$collection_name[which(tax.dat_Apt$lagerstatten=="conservation")])
length(L_names_Apt)

L_names_Cam <- unique(tax.dat_Cam$collection_name[which(tax.dat_Cam$lagerstatten=="conservation")])
length(L_names_Cam)

L_names_Toa <- unique(tax.dat_Toa$collection_name[which(tax.dat_Toa$lagerstatten=="conservation")])
length(L_names_Toa)


# 2. Collections per latitude ------------------------------------------------

## Yesterday, you had a look with alpha diversity ('local richness') 
##    When visualised, alpha diversity can also provide more insights into sampling patterns
##    especially as it adds a spatial element as opposed to just temporal patterns.
## Let's create a plot to see where the collections across time and palaeolatitude.
##    We'll also colour our plot according to the number of taxa in each collection

## There is evidence to suggest that alpha diversity is not as strongly affected by sampling biases
##    as gamma (or 'global') diversity. For a more sophisticated way to calculate alpha diversity by
##    treating taxonomically indeterminate occurrences as valid, see the method described in 
##    Close et al. (2019) - code available here: https://github.com/emmadunne/local_richness

## Let's get our data set up:
lat_data <- occurrences_sp # rename object to keep the original separate

## Create new column for mid_ma
lat_data$mid_ma <- (lat_data$max_ma + lat_data$min_ma)/2 

## Next, we'll need to count the number of taxa per collection (i.e. their frequency):
taxa_freqs <- count(lat_data, collection_no)

## Subset lat_data to only the columns we need:
lat_data <- lat_data %>% 
  select(collection_no, paleolat, paleolng, mid_ma) %>% 
  distinct() %>% na.omit()

## And add the frequency information:
lat_data <- left_join(taxa_freqs, lat_data, by = "collection_no")

## Before we plot, let's order the frequencies and remove any NAs that have crept in:
lat_data <- lat_data %>% arrange(n) %>% na.omit()

## Take a look:
View(lat_data)


## Set up our ggplot layers
lat_plot <- ggplot(data = lat_data, aes(x = mid_ma, y = paleolat, colour = n)) +
  geom_vline(xintercept = int_boundaries, lty = 2, col = "grey90") +
  geom_hline(yintercept = 0, colour = "grey10") +
  scale_color_viridis(trans = "log", breaks = c(1, 2, 6, 14), direction = -1, option = "D") + # set the break= to match your richness data
  #scale_y_continuous(labels = function(x) format(x, width = 5), limits = c(-70, 70), breaks = seq(from = -60, to = 60, by = 20)) +
  scale_x_reverse(breaks = int_boundaries) + 
  theme_minimal() + 
  theme(legend.direction = "vertical", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Palaeolatitude (º)") +
  geom_point(size = 4, alpha = 0.5) # (alpha sets point transparency)
lat_plot # call to plot window


## Set dimensions and save plot (as pdf)
ggsave(plot = lat_plot,
       width = 20, height = 10, dpi = 500, units = "cm", 
       filename = "./plots/lat_alpha_div.pdf", useDingbats=FALSE)

# find data point with exceptional number and low lat
# minlat <- min(abs(lat_data$paleolat))
# mincol <- lat_data[which(abs(lat_data$paleolat)== minlat),]
# 
# # exclude that point
# latdat <- lat_data[-which(abs(lat_data$paleolat)== minlat),]
# lat_data[which(abs(latdat$paleolat)== min(abs(latdat$paleolat))),]

lat_dat_ord <- lat_data[order(abs(lat_data$paleolat),decreasing = F),]
loc <- lat_dat_ord [6,]

# # location today
location <- select(occurrences_sp, collection_no, collection_name, lat, lng,
                   paleolat, paleolng, max_ma, min_ma)
location[which(location$paleolat == loc$paleolat),] # Brazil


# 3. Regression -----------------------------------------------------------

## We can also apply some simple statistics to better quantify the
##    relationship between raw richness and sampling proxies/effort
## Let's do some simple regression plots:

## Raw richness vs. collections
reg_colls <- ggplot(proxy_counts, aes(x=count_taxa, y=count_colls)) + 
  geom_point(shape=17, size = 6, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1")  +
  theme_minimal()
reg_colls

## Raw richness vs. formations
reg_forms <- ggplot(proxy_counts, aes(x=count_taxa, y=count_formations)) + 
  geom_point(shape=16, size = 5, colour = "#30C430")+
  geom_smooth(method=lm, colour = "#0A6B09", fill = "#B0ECB0") +
  theme_minimal()
reg_forms


## Let's quantify these relationships through a linear model:

## Raw richness vs. collections
lm_colls = lm(count_colls ~ count_taxa, proxy_counts)
summary(lm_colls) # summary of results

lm_forms = lm(count_formations ~ count_taxa, proxy_counts)
summary(lm_forms)

## We can test the sensitivity of the data to well-sampled intervals
## In this case, we could remove the Norian and/or the Carnian and re-run 
##    the analyses using the code above

## Take out row 6, which is the Norian
proxy_counts_N <- proxy_counts[-6,]
## Take out row 6 + 7 (Norian + Carnian)
proxy_counts_NC <- proxy_counts[-c(6,7),]

## Are there any changes to your results? 



# 4. Collector's curve -------------------------------------------------------

## Collector's curves, also known as a species accumulation curves, can also tell us about 
##    sampling: they display the cumulative number of taxa as a function of the cumulative 
##    effort (i.e. sampling proxy)

## First, get our data into shape:
## Table the number of each species per collection (abundance)
abun_data <- table(occurrences_sp$collection_no, occurrences_sp$accepted_name)
## If we want a presence/absence table, we can add this step:
#abun_data[abun_data > 0] <- 1

## Turn this into a matrix:
abun_matrix <- matrix(abun_data, ncol = length(unique(occurrences_sp$accepted_name)))
colnames(abun_matrix) <- colnames(abun_data) # add the column names back in for when we need to check anything
rownames(abun_matrix) <- rownames(abun_data) # same for the row names

## Using the vegan package, we can make the collectors or species accumulation curve
## Check out the help file for specaccum() to find out more about the methods
sp_accum <- specaccum(abun_matrix, method = "collector")

## Plot the curve in base R - you can make this pretty in ggplot if you prefer!
plot(sp_accum, ci.type = "poly", col = "#0E6A8A", lwd = 2, ci.lty = 0, ci.col = "#5CBCDD")

