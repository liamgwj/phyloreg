# LJ 2022-03-07 prep inputs and source model fitting script
## to be used interactively!

library(ape)

setwd("/home/liam/Documents/MSc/analysis/phyloreg")

# choose dataset ---------------------------------------------------------------

# Parker
topdir <- "datasets/Parker/"

dbfile <- "cleaned_Parker_database_2022-02-28_p067"

db <- read.csv(paste0(topdir,"cleaned_data/database_subsets/",
                      dbfile, ".csv"),
               row.names = 1)
phy <- read.tree(paste0(topdir, "cleaned_data/", "QJ_Parker_2022-03-09.nwk"))

# Robles

# simulated
topdir <- "datasets/sim/"

dbfile <- "sim-20220315T121340_db_p067"

db <- read.csv(paste0(topdir,"simulated_data/database_subsets/",
                      dbfile, ".csv"),
               row.names = 1)

phy <- read.tree(paste0(topdir,"simulated_data/","sim-20220315T121340_phy.nwk"))


# if necessary, remove columns/tips to ensure perfect intersection of db and phy

db <- db[,which(colnames(db)%in%phy$tip.label)]

phy <- keep.tip(phy, colnames(db))

#if(min(colSums(db))==0){
#db <- db[which(rowSums(db)==0) ,]
#}


# source model script ----------------------------------------------------------

source("/home/liam/Documents/MSc/analysis/phyloreg/scripts/logreg_generic.r")

