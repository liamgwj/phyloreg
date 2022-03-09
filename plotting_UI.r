# LJ 2022-03-07 set up for plotting and source plotting script
## to be used interactively!

# choose dataset/model coefficient files to plot -------------------------------

# Parker 2015 -------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

coef_1 <- read.csv("output/model_coefficients/parker-QJ-all_coef-av_2022-03-09T12-31-19.csv")
coef_2 <- read.csv("output/model_coefficients/parker-QJ-half_coef-av_2022-03-09T12-33-37.csv")


# Robles 2017 -------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Robles")

coef_1 <- read.csv("output/model_coefficients/robles-all_coef-av_2022-03-08T09-57-30.csv")
coef_2 <- read.csv("output/model_coefficients/robles-half_coef-av_2022-03-08T09-58-43.csv")


# format list of coefficients --------------------------------------------------

coef_lst <- vector("list", length = nrow(coef_2))

names(coef_lst) <- coef_2$Pest

for(i in 1:nrow(coef_2)){
coef_lst[[i]] <- rbind(coef_1[which(coef_1$Pest == coef_2$Pest[i]) ,],
                       coef_2[i,])
}


# source plotting script -------------------------------------------------------

# set ID for file naming
ID <- "parker-QJ-ah"

source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")

