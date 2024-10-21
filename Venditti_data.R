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
# 
# ******************************************************

# load in data
eff.dat <- read.table("Data/Input/Venditti_avgData.txt")

# new header
names(eff.dat)
names(eff.dat) <- eff.dat[1,]
names(eff.dat)

eff.dat <- eff.dat[-1,]

write.csv2(eff.dat,"Data/Output/Species_eff_data.csv", row.names = FALSE)
