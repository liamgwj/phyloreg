# LJ 2022-02-28 generic model prediction and plot generation script

# requires objects:
# 'coef', as produced by 'logreg_generic.r'
# 'dataSource' (to be used when naming plot files)
# to be present in the environment

# average coefficients for each pest
coef_av <- aggregate(coef[,-which(colnames(coef)%in%c("Pest", "FocalHost"))],
                     list(Pest = coef$Pest),
                     mean)

# generate per-pest predictions for plotting curves

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

