# LJ 2022-02-18 clean up pest/host database and phylogeny from Parker et al.
# 2015, including addition of many missing genera to phylogeny

# SETUP #######################################################################

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

library(ape)
library(phytools)
library(phangorn)


# DATA IMPORT #################################################################

# anonymized pest database from Gilbert 2012/Parker 2015
db <- read.csv("indata/eva265-sup-0001-datas1.csv")

# host phylogeny from Parker
phy <- read.tree("indata/41586_2015_BFnature14372_MOESM285_ESM.txt")

# plant family-genus table according to TPL
tplfg <- read.csv("indata/tpl_fg.csv")


# DATA CLEANUP ################################################################

# pest database -----------------------------------------------------

# there is one NA value at [269,214] - not sure why
# set to 0 for now
db[269,214] <- 0

# some pests (rows) have no recorded hosts (1 values - 0's only)
# take out zero-host pest rows
db <- db[-which(rowSums(db[,2:ncol(db)])==0),]

# get list of included genera
dbgen <- colnames(db)[-which(colnames(db)=="Pestcode")]

# plant family-genus table -------------------------------------------

# match format to db
tplfg$family <- tolower(tplfg$Family)
tplfg$genus <- tolower(tplfg$Genus)

# subset to host genera
dbfg <- tplfg[tplfg$genus%in%dbgen,]

# correct for taxonomic differences between TPL and APGIII
dbfg$fam_cor <- dbfg$family
dbfg$fam_cor[which(dbfg$family == "compositae")] <- "asteraceae"
dbfg$fam_cor[which(dbfg$family == "leguminosae")] <- "fabaceae"

# some genus names are shared between families
# we assume that the genus name in the pest db refers to one unique pairing
# need to figure out which pairings are correct
# Ditrichaceae are mosses that share a genus name with Poaceae
# confident about removing Ditrichaceae/Cynodon, others need checking

# for now, arbitrarily keeping (mostly) first case
dbfg <- dbfg[-c(which(dbfg$Family=="Ditrichaceae" & dbfg$Genus=="Cynodon"),
                which(dbfg$Family=="Clusiaceae" & dbfg$Genus=="Calophyllum"),
                which(dbfg$Family=="Moraceae" & dbfg$Genus=="Urtica"), #not 1st
                which(dbfg$Family=="Rhamnaceae" & dbfg$Genus=="Ziziphus")
                ),]

# phylogeny -----------------------
# make ultrametric
phy <- force.ultrametric(phy, method = "extend")
## this is a crappy workaround - need to properly trace problem

# IDENTIFY GENERA TO ADD ######################################################

# length(dbfg$genus) # 213 host genera in database
# length(intersect(dbfg$genus, phy$tip.label)) # 78 genera in phy as tips

# identify the 135 genera to add to phy
nottips <- dbfg[-which(dbfg$genus%in%phy$tip.label),]

# of these, 132 have their family in phy as a tip
famtips <- nottips[which(nottips$fam_cor%in%phy$tip.label),]

# id the remaining 3 genera
nofamtips <- nottips[-which(nottips$fam_cor%in%phy$tip.label),]

# tally the number of genera per family in database
genperfam <- vector(length=length(unique(dbfg$fam_cor)))
names(genperfam) <- unique(dbfg$fam_cor)

for(i in seq_along(genperfam)){
tmp <- dbfg[which(dbfg$fam_cor == names(genperfam)[i]),]
genperfam[i] <- length(tmp$gen)
}

# id single-genus families
sgfam <- names(genperfam[which(genperfam == 1)])

# divide genera with families as tips by whether or not they are singletons
sfamtips <- famtips[which(famtips$fam_cor %in% sgfam),]
mfamtips <- famtips[-which(famtips$fam_cor %in% sgfam),]

# addition scenarios
# CASE 1: family is phy as tip, single genus in family (sfamtips)
# change family tip name to genus name

# CASE 2: family in phy as tip, multiple genera (mfamtips)
# add all genera as a polytomy at 1/2 way point of edge leading to family tip

# CASE 3: family not in phy as tip (nofamtips)
# will deal with these on a case-by-case basis (only 3 genera affected)


# PHYLOGENETIC ADDITIONS ######################################################

# CASE 1 ---------------------------------------------------

# swap family tip names for genus names
for(i in seq_along(sfamtips$fam_cor)){
phy$tip.label[which(phy$tip.label==sfamtips$fam_cor[i])] <- sfamtips$genus[i]
}

# CASE 2 ---------------------------------------------------

# get lengths of all terminal edges of phy (from phytools blog)
edgelengths <- setNames(phy$edge.length[sapply(1:Ntip(phy),
                                               function(x,y){which(y==x)},
                                               y=phy$edge[,2])],
                        phy$tip.label)

# specify "addition length" for each genus (1/2 family's terminal edge length)
addlengths <- edgelengths[which(names(edgelengths) %in%
                                                    unique(mfamtips$fam_cor))]
addlengths <- addlengths / 2

# match ordering to family list
addlengths <- addlengths[match(unique(mfamtips$fam_cor), names(addlengths))]

# generate polytomous subtrees
fam_subtrees <- vector(length = length(unique(mfamtips$fam_cor)))
names(fam_subtrees) <- unique(mfamtips$fam_cor)

for(i in 1:length(fam_subtrees)){
tmp <- mfamtips[which(mfamtips$fam_cor == unique(mfamtips$fam_cor)[i]),]
fam_subtrees[i] <- paste0("(")
for(j in 1:length(tmp$genus)){
if(j == length(tmp$genus)){
fam_subtrees[i] <- paste0(fam_subtrees[i],
                          tmp$genus[j],
                          ":",
                          addlengths[i])
}else{
fam_subtrees[i] <- paste0(fam_subtrees[i],
                          tmp$genus[j],
                          ":",
                          addlengths[i],
                          ",")
}}
fam_subtrees[i] <- paste0(fam_subtrees[i],
                          "):",
                          addlengths[i],
                          ";")
}

# bind to main tree
for(i in 1:length(fam_subtrees)){
phy <- bind.tree(phy,
                 read.tree(text=fam_subtrees[i]),
                 phy$edge[match(match(names(fam_subtrees)[i],
                                      phy$tip),
                                phy$edge[,2]),
                          1])
}



####
###
# plot(extract.clade(phy, getMRCA(phy,c("rosaceae","malus")))) 
##
###

# CASE 3 ---------------------------------------------------

# 3 genera from 2 families to be handled case-by-case

# genus quercus is present in phy, but represented by many species-level tips
# change name of first species to "quercus"
phy$tip.label[grep("quercus", phy$tip.label)[1]] <- "quercus"

# family solanaceae is present as a labeled node, but neither genus is present
# database includes other genera from family that are already present in phy
# add genera as polytomy at base of family
sl <- edgelengths[which(names(edgelengths)=="schizanthus")]

phy <- bind.tip(phy,
                "lycopersicon",
                sl,
                which(phy$node.label=="solanaceae") + Ntip(phy))

phy <- bind.tip(phy,
                "solanum",
                sl,
                which(phy$node.label=="solanaceae") + Ntip(phy))


# WRITE OUTPUT ################################################################

# Final misc rearranging - intergrate this earlier
# db
# set row names to pest names
rownames(db) <- db$Pestcode
# remove pest name column
db <-  db[, !names(db)%in%"Pestcode"]
# remove pests that infect only one host
db <- db[which(rowSums(db)>1),]
# phy
# drop tips not in pest/host database
phy <- keep.tip(phy, colnames(db))


write.tree(phy,
           paste0("output/cleaned_Parker_phylogeny_",
                  Sys.Date(),
                  ".nwk"))

write.csv(db,
          paste0("output/cleaned_Parker_database_",
                 Sys.Date(),
                 ".csv"),
          row.names = TRUE)

