# LJ 2022-03-14 sample a subset of the host/pest matrix

# env holds: db

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

if(!dir.exists(paste0(dir, "/database_subsets"))){
dir.create(paste0(dir, "/database_subsets"))
}

write.csv(db,
          paste0(dir, "/database_subsets/", file,
                 "_p", gsub("\\.", "", round(prop, 2)), ".csv"))

