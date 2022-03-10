# LJ 2022-03-07 prep inputs and source model fitting script
## to be used interactively!

library(ape)

# choose dataset ---------------------------------------------------------------

# Parker 2015 --------------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

db_safe <- read.csv("output/cleaned_Parker_database_2022-02-28.csv",
                    row.names = 1)

#phy_safe <- read.tree("output/cleaned_Parker_phylogeny_2022-02-28.nwk")

phy_safe <- read.tree("output/QJ_Parker_2022-03-09.nwk")

db_safe <- db_safe[,colnames(db_safe)%in%phy_safe$tip.label]

# Robles 2017 --------------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Robles")

db_safe <- read.csv("output/cleaned_Robles_database_2022-03-07.csv",
                    row.names = 1)

phy_safe <- read.tree("indata/host.tree")


# subset as desired ------------------------------------------------------------

# all pests ----------------------------

db <- db_safe
phy <- phy_safe

# random subset of pests ---------------
# choose proportion of pests to keep
prop <- 1/2

# remove 1-prop rows from db
db <- db_safe[sample(1:nrow(db_safe), prop*nrow(db_safe)),]

# if any hosts have had all their pests removed, remove them
if(min(colSums(db)) == 0){
db <- db[, -which(colSums(db) == 0)]
}

# inverse
db <- db_safe[,sample(1:ncol(db_safe), prop*ncol(db_safe))]

db <- db[-which(rowSums(db)==0),]


# prune phy to remaining hosts
phy <- keep.tip(phy_safe, colnames(db))

# single pest type ---------------------
# Parker types: Bacterium Fungus Insect Mite Mollusk Nematode Phytoplasma
# Viroid Virus Weed

db <- db_safe[grep("Weed", rownames(db_safe)),]

if(min(colSums(db)) == 0){
db <- db[, -which(colSums(db) == 0)]
}

phy <- keep.tip(phy_safe, colnames(db))


# source model script ----------------------------------------------------------

# set ID for file naming
ID <- "parker-QJ-halfhosts"

source("/home/liam/Documents/MSc/analysis/phyloreg/logreg_generic.r")

