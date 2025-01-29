########### script for filtering PBDB quickly



## Download names and ages of time intervals from the PBDB:
intervals_all <- read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")
View (intervals_all) # take a look

## For the rest of this session, we're going to focus on the Late Triassic-Early Jurassic interval
## Make a vector of stage names that we are interested in:
interval_names <- c("Carnian", "Norian", "Rhaetian", # Late Triassic
                    "Hettangian", "Sinemurian", "Pliensbachian", "Toarcian") # Early Jurassic

## Select these intervals from the full PBDB intervals dataset:
intervals <- filter(intervals_all, interval_name %in% interval_names)

## Pare this down to just the 3 columns we'll need:
intervals <- select(intervals, interval_name, early_age, late_age)

## For ease of use later, let's rename the age columns to match the occurrence data:
intervals <- rename(intervals, "max_ma" = "early_age", "min_ma" = "late_age")

## And finally, calculate the midpoint for each interval and add it to a new (4th) column
intervals$mid_ma <- (intervals$min_ma + intervals$max_ma)/2

## Take a peep:
View(intervals) # open as new tab

## Save a copy as a .csv file - Note: your file path will differ!
write_csv(intervals, "./data/intervals_Car_Tor.csv")