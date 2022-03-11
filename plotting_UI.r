# LJ 2022-03-07 set up for plotting and source plotting script
## to be used interactively!

library(ape)

# choose dataset/model coefficient files to plot -------------------------------

# Parker 2015 -------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

coef_1 <- read.csv("output/model_coefficients/parker-QJ-all_coef-av_2022-03-09T12-31-19.csv")
coef_2 <- read.csv("output/model_coefficients/parker-QJ-halfhosts_coef-av_2022-03-10T15-38-24.csv")

# for phylo
coef_2a <- read.csv("output/model_coefficients/parker-QJ-halfhosts_coef_2022-03-10T15-38-24.csv")
phy <- read.tree("output/QJ_Parker_2022-03-09.nwk")
#phy <- read.tree("output/cleaned_Parker_phylogeny_2022-02-28.nwk")
db <- read.csv("output/cleaned_Parker_database_2022-02-28.csv",
               row.names = 1)

# Robles 2017 -------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Robles")

coef_1 <- read.csv("output/model_coefficients/robles-all_coef-av_2022-03-08T09-57-30.csv")
coef_2 <- read.csv("output/model_coefficients/robles-half_coef-av_2022-03-08T09-58-43.csv")

coef_2a <- read.csv("output/model_coefficients/robles-half_coef_2022-03-08T09-58-43.csv")

phy <- read.tree("indata/host.tree")

db <- read.csv("output/cleaned_Robles_database_2022-03-07.csv",
               row.names = 1)


# format list of coefficients --------------------------------------------------

coef_lst <- vector("list", length = nrow(coef_2))

names(coef_lst) <- coef_2$Pest

for(i in 1:nrow(coef_2)){
coef_lst[[i]] <- rbind(coef_1[which(coef_1$Pest == coef_2$Pest[i]) ,],
                       coef_2[i,])
}

# phylogeny plotting ----------------------------------------------------------




# curve plotting ---------------------------------------------------------------

# set ID for file naming
ID <- "parker-QJ-ah"

source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")

