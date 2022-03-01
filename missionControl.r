# LJ 2022-02-27 Central script for setting up and sourcing the various modules
# this script is meant to be used interactively

# SETUP ########################################################################

library(ape)


# PARKER 2015 ##################################################################

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

# compatibility matrix
db_safe <- read.csv("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/cleaned_Parker_database_2022-02-28.csv", row.names = 1)

# phylogeny
phy_safe <- read.tree("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/cleaned_Parker_phylogeny_2022-02-28.nwk")

# all pests ------------------------------------------------------------

db <- db_safe
phy <- phy_safe

# source logistic regression model script
source("/home/liam/Documents/MSc/analysis/phyloreg/logreg_generic.r")

# read in model output coefficients
coef_safe <- read.csv("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/All/coefficients_2022-02-28.csv")

# identifier for plot naming
dataSource <- "Parker-all"

# source model prediction/plotting script
source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")

# get pest categories
pestTypes <- unique(gsub("[0-9]","", coef$Pest))

# loop to make by-catergory plots
for(k in seq_along(pestTypes)){
# name datasource
dataSource <- paste0("Parker-", pestTypes[k])
# subset coefficients
coef <- coef_safe[grep(pestTypes[k], coef_safe$Pest) ,]
# source plotting script
source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")
}

# pests by type -------------------------------------------------------

# get pest categories
pestTypes <- unique(gsub("[0-9]", "", rownames(db_safe)))

# models
# for each pest type:
for(k in seq_along(pestTypes)){
# subset db to only one pest type
db <- db_safe[grep(pestTypes[k], rownames(db_safe)), ]
# remove hosts without any associations
db <- db[, which(colSums(db)>0)]
# prune phy to only remaining hosts
phy <- keep.tip(phy_safe, colnames(db))
# create and move into output directory
dir.create(paste0("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/", pestTypes[k]))
setwd(paste0("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/", pestTypes[k]))
# source logistic regression model script
source("/home/liam/Documents/MSc/analysis/phyloreg/logreg_generic.r")
}

# plotting
for(k in seq_along(pestTypes)){
# navigate
setwd(paste0("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/", pestTypes[k]))
# name datasource
dataSource <- paste0("Parker-", pestTypes[k])
# subset coefficients
coef <- read.csv("coefficients_2022-02-28.csv")
# source plotting script
source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")
}

