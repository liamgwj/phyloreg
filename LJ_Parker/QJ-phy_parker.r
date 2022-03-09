# LJ 2022-03-08 try parker analysis with QJ phylogeny

library(ape)

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")


phy_safe <- read.tree("indata/PhytoPhylo.tre")

db <- read.csv("output/cleaned_Parker_database_2022-02-28.csv",
                row.names = 1)

# 211 of 213 genus names from db are present in phy - awesome!
# however, we will still have to account for duplicated names across families

phy <- keep.tip(phy_safe, which(tolower(gsub("_.*", "", phy_safe$tip.label)) %in% colnames(db)))

tax <- data.frame(tip=phy$tip.label)

tax$Genus <- gsub("_.*", "", tax$tip)

# plant family-genus table according to TPL
tplfg <- read.csv("indata/tpl_fg.csv")

tax <- merge(tax, tplfg)

tax <- tax[-c(which(tax$Family=="Ditrichaceae" & tax$Genus=="Cynodon"),
                which(tax$Family=="Clusiaceae" & tax$Genus=="Calophyllum"),
                which(tax$Family=="Moraceae" & tax$Genus=="Urtica"),
                which(tax$Family=="Rhamnaceae" & tax$Genus=="Ziziphus")
                ),]

# collapse genera to single tips
# identify first tip-level child of each genus MRCA

phy <- keep.tip(phy, tax$tip)


gentips <- vector("list", length = length(unique(tax$Genus)))

names(gentips) <- unique(tax$Genus)

mrcas <- vector(length = length(unique(tax$Genus)))

keeptips <- vector(length = length(unique(tax$Genus)))


for(i in seq_along(unique(tax$Genus))){

gentips[[i]] <- tax$tip[which(tax$Genus == unique(tax$Genus)[i])]


if(length(gentips[[i]]) == 1){

keeptips[i] <- gentips[[i]]

mrcas[i] <- NA

}else{

mrcas[i] <- getMRCA(phy, gentips[[i]])


if(!is.na(mrcas[i])){

keeptips[i] <- extract.clade(phy, mrcas[i])$tip.label[1]

}}}


phy <- keep.tip(phy, keeptips)

phy$tip.label <- unique(tax$Genus)










