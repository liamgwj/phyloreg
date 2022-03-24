# LJ 2022-03-07 isolate Robles cleaning from messy script in prep for generic workflow

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Robles")

library(ape)

# host phylogenetic tree built with phylomatic
phy <- read.tree("host.tree")

# full beetle/host tree compatibility database
db_full <- read.csv("full_taxa_beetle_plant.csv", header = T)

## note that db_full includes associations to species, but phy includes
## tips only at genus level


# host/pest database ----------------------------------------------

# convert paired columns to binary matrix
db <- table(db_full[,c(7,13)])

# remove values >1
db <- ifelse(db > 1, 1, db)

# this database looks good - all hosts have at least one pest, all pests have at least two hosts

# however, 8 host genera are not represented in phy - remove them from db
db <- db[,-which(!colnames(db)%in%phy$tip.label)]

# write to file
write.csv(db,
          paste0("cleaned_Robles_database_",
                 Sys.Date(),
                 ".csv"),
          row.names = TRUE)

