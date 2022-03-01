# modified from
# LJ 2022-02-28 generic model prediction and plot generation script


pestTypes <- c("Bacterium", "Fungus", "Insect", "Mite", "Mollusk", "Nematode",
               "Phytoplasma", "Viroid", "Virus", "Weed")


# read in whole dataset model coefs
coef_all <- read.csv("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/All/coefficients_2022-02-28.csv")

# average coefficients across runs for each pest species
coef_av <- aggregate(coef_all[,-which(colnames(coef_all)%in%c("Pest", "FocalHost"))],
                     list(Pest = coef_all$Pest),
                     mean)


for(i in seq_along(pestTypes)){

coef_it <- read.csv(paste0("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/", pestTypes[i], "/coefficients_2022-02-28.csv"))

# average coefficients for each pest
coef_it_av <- aggregate(coef_it[,-which(colnames(coef_it)%in%c("Pest", "FocalHost"))],
                     list(Pest = coef_it$Pest),
                     mean)

# list to hold predictions
curve_data_it <- vector("list", length = length(coef_it_av$Pest))
names(curve_data_it) <- coef_it_av$Pest

# loop to fill list
for(j in seq_along(coef_it_av$Pest)){
# generate predicted values at a range of PDs
logit_it <- coef_it[i, "intercept"] + coef_it[i, "slope"] * log10(seq(0, 650, 1) + 1)
# convert to probabilities
prob_logit_it <- (exp(logit_it)/(1 + exp(logit_it)))
# if predicted value was too big, Inf division has given NA - convert these
# values to 1 and store
curve_data_it[[j]] <- ifelse(is.na(prob_logit_it), 1, prob_logit_it)
}

# all data predictions ----------------------------

all_it <- coef_av[grep(pestTypes[i], coef_av$Pest), ]

# list to hold predictions
curve_data <- vector("list", length = length(all_it$Pest))
names(curve_data) <- all_it$Pest

# loop to fill list
for(j in seq_along(all_it$Pest)){
# generate predicted values at a range of PDs
logit <- all_it[j, "intercept"] + all_it[j, "slope"] * log10(seq(0, 650, 1) + 1)
# convert to probabilities
prob_logit <- (exp(logit)/(1 + exp(logit)))
# if predicted value was too big, Inf division has given NA - convert these
# values to 1 and store
curve_data[[j]] <- ifelse(is.na(prob_logit), 1, prob_logit)
}

# open plotting connection
jpeg(paste0(pestTypes[i], "_predSuit.jpg"),
     width = 85, height = 85, units = 'mm', res = 300)

# base plot to provide axis labels, actual lines are omitted
matplot(curve_data[[1]], type="n",
        ylab="Probability of sharing pest species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F)

# mean line (all)
lines(colMeans(data.frame(t(sapply(curve_data, c)))), col="firebrick")

# mean line (it)
lines(colMeans(data.frame(t(sapply(curve_data_it, c)))), col="blue")

# x axis (below)
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis = "black", las = 1, cex.axis = 0.5)
# y axis (left)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis = "black", las = 1, cex.axis = 0.5)

# legend
legend(300, 1,
       legend = c("all pests",
                  pestTypes[i]),
       lty = c(1, 1),
       col = c("firebrick", "blue"),
       lwd = c(3, 3),
       cex = 0.3,
       bty = "n")

# close connection
dev.off()

}


