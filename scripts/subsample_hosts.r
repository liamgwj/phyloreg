# LJ 2022-03-14 sample a subset of the host/pest matrix

# read in data to subsample ---------------------------------------------------

dir <- "indata/simulated/"
runID <- "sim-20220328T124900"

db <- read.csv(paste0(dir, runID, "_db.csv"), row.names = 1)


# subset of hosts -------------------------------------------------------------

# set sampling proportion
prop <- 2/3


for(i in seq_along(rownames(db))){

# for a given pest
focal <- rownames(db)[i]

# get all potential hosts
pothosts <- colnames(db)

# remove known hosts
hostind <- which(db[focal,]==1)
pothosts <- pothosts[-hostind]

# remove 1-prop potential hosts
pothosts <- sample(pothosts, prop*length(pothosts), replace = FALSE)

# re-add known hosts
keepers <- c(pothosts, colnames(db)[hostind])

# set db values for removed potential hosts to NA
db[focal,] <- ifelse(!colnames(db)%in%keepers, NA, db[focal,])

}

# write out -------------------------------------------------------------------

write.csv(db,
          paste0("output/", runID,
                 "_p", gsub("\\.", "", round(prop, 2)), "_db.csv"))

