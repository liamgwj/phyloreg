# LJ 2022-02-18 get tpl family-genus table

# reading Parker, the missing genera were added as polytomies at the family
# node... but there's no taxonomy file included! just a statement that they
# "used TPL"... ugh

library(taxize)

# download the whole plant list to a folder called 'TPL' on the desktop
# one CSV file per family
tpl_get("~/Desktop/TPL")

# run bash script to extract family/genus pairs

# read in result
tplfg <- read.csv("~/Desktop/tpl_famgen.csv")

# unique cases
tplfg <- unique(tplfg)

# save
write.csv(tplfg, "tpl_fg.csv", row.names=FALSE)
