# modified from
# LJ 2022-02-28 generic model prediction and plot generation script

# read in model coefs
coef_all <- read.csv("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker/output/All/coefficients_2022-02-28.csv")

# average coefficients for each pest
coef_av <- aggregate(coef[,-which(colnames(coef)%in%c("Pest", "FocalHost"))],
                     list(Pest = coef$Pest),
                     mean)

# list to hold predictions
curve_data <- vector("list", length = length(coef_av$Pest))
names(curve_data) <- coef_av$Pest

# loop to fill list
for(i in seq_along(coef_av$Pest)){
# generate predicted values at a range of PDs
logit <- coef[i, "intercept"] + coef[i, "slope"] * log10(seq(0, 650, 1) + 1)
# convert to probabilities
prob_logit <- (exp(logit)/(1 + exp(logit)))
# if predicted value was too big, Inf division has given NA - convert these
# values to 1 and store
curve_data[[i]] <- ifelse(is.na(prob_logit), 1, prob_logit)
}










# open plotting connection
jpeg(paste0(dataSource, "_predSuit.jpg"),
     width = 85, height = 85, units = 'mm', res = 300)

# base plot to provide axis labels, actual lines are omitted
matplot(curve_data[[1]], type="n",
        ylab="Probability of sharing pest species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F)

# add one line per model run
lapply(curve_data, function(x){lines(x, col="grey82")})
# mean line
lines(colMeans(data.frame(t(sapply(curve_data, c)))), col="firebrick")

# x axis (below)
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis = "black", las = 1, cex.axis = 0.5)
# y axis (left)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis = "black", las = 1, cex.axis = 0.5)

# close connection
dev.off()





pestTypes <- c("Bacterium", "Fungus", "Insect", "Mite", "Mollusk", "Nematode",
               "Phytoplasma", "Viroid", "Virus", "Weed")

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
curve_data_it[[i]] <- ifelse(is.na(prob_logit_it), 1, prob_logit_it)
}


}



