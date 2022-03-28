# LJ 2022-03-07 set up for plotting and source plotting script
## to be used interactively!

library(ape)

#setwd("/home/liam/Documents/MSc/analysis/phyloreg")

# choose dataset/model coefficient files to plot -------------------------------

# parker
#topdir <- "datasets/Parker/"
#
#coef1id <- "cleaned_Parker_database_2022-02-28_GroupCoefs"
#
#coef2id <- "cleaned_Parker_database_2022-02-28_p067_GroupCoefs"
#
#coef_1 <- read.csv(paste0(topdir, "model_coefficients/", coef1id, ".csv"))
#
#coef_2 <- read.csv(paste0(topdir, "model_coefficients/", coef2id, ".csv"))
#
#phy <- read.tree(paste0(topdir, "cleaned_data/",
#                        "QJ_Parker_2022-03-09.nwk"))
#
#db <- read.csv(paste0(topdir, "cleaned_data/",
#                      "cleaned_Parker_database_2022-02-28.csv"),
#               row.names = 1)
#
#db_2 <- read.csv(paste0(topdir, "cleaned_data/database_subsets/",
#                        "cleaned_Parker_database_2022-02-28_p067.csv"),
#                 row.names = 1)
#
#db <- db[,which(colnames(db)%in%phy$tip.label)]
#phy <- keep.tip(phy, colnames(db))
#if(min(colSums(db))==0){
#db <- db[which(rowSums(db)==0) ,]}


# robles



# simulated
runID <- "sim-20220328T124900"

subID <- "_p067"

coef_1 <- read.csv(paste0("output/", runID, "_coef-av.csv"))

coef_2 <- read.csv(paste0("output/", runID, subID, "_coef-av.csv"))

phy <- read.tree(paste0("indata/simulated/", runID, "_phy.nwk"))

db <- read.csv(paste0("indata/simulated/", runID, "_db.csv"),
               row.names = 1)

db_2 <- read.csv(paste0("indata/simulated/", runID, subID, "_db.csv"),
                 row.names = 1)


# format list of coefficients --------------------------------------------------

coef_lst <- vector("list", length = nrow(coef_2))

names(coef_lst) <- coef_2$Pest

for(i in 1:nrow(coef_2)){
coef_lst[[i]] <- rbind(coef_1[which(coef_1$Pest == coef_2$Pest[i]) ,],
                       coef_2[i,])
}

# phylogeny plotting ----------------------------------------------------------

source("/home/liam/Documents/MSc/analysis/phyloreg/r/scripts/functions/plot_phylogeny.r")


# curve plotting ---------------------------------------------------------------

source("/home/liam/Documents/MSc/analysis/phyloreg/r/scripts/functions/plot_curves.r")

