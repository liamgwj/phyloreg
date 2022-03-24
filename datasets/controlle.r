
setwd("/home/liam/Documents/MSc/analysis/phyloreg/datasets")

# read in data to subsample ---------------------------------------------------

# Parker 2015
dir <- "Parker/cleaned_data/"
file <- "cleaned_Parker_database_2022-02-28"

# simulated
#dir <- "sim/simulated_data/"
#file <- "sim-20220315T121340_db"


db <- read.csv(paste0(dir, file, ".csv"), row.names = 1)


