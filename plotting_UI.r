# LJ 2022-03-07 set up for plotting and source plotting script
## to be used interactively!

# choose dataset/model coefficient files to plot -------------------------------

# Robles 2017 -------------------

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Robles")

coef_1 <- read.csv("output/model_coefficients/coef-av_2022-03-07T14-56-18.csv")
coef_2 <- read.csv("output/model_coefficients/coef-av_2022-03-07T15-00-21.csv")


# format list of coefficients --------------------------------------------------

coef_lst <- vector("list", length = nrow(coef_2))

names(coef_lst) <- coef_2$Pest

for(i in 1:nrow(coef_2)){
coef_lst[[i]] <- rbind(coef_1[which(coef_1$Pest == coef_2$Pest[i]) ,],
                       coef_2[i,])
}


# source plotting script -------------------------------------------------------

source("/home/liam/Documents/MSc/analysis/phyloreg/predict_plot_generic.r")

