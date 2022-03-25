# LJ 2022-02-28 generic model prediction and plot generation script

# generate one plot per pest, each with one line per model

# 'coef_lst' is a list with one entry per pest, each entry containing one row per model run to compare

# dirname


# generate per-pest predictions for plotting curves

logit_lst <- vector("list", length = length(coef_lst))
names(logit_lst) <- names(coef_lst)

for(i in 1:length(coef_lst)){

for(j in 1:nrow(coef_lst[[i]])){

if(j == 1){

logit_lst[[i]] <- list(coef_lst[[i]][j,"intercept"] +
                       (coef_lst[[i]][j,"slope"] * log10(seq(0,650,1) + 1)))
}else{

logit_lst[[i]][j] <- list(coef_lst[[i]][j,"intercept"] +
                          (coef_lst[[i]][j,"slope"] * log10(seq(0,650,1) + 1)))
}}}

# convert to predicted probability
prob_lst <- rapply(logit_lst, function(x){(exp(x)/(1 + exp(x)))}, how='list')

# if predicted value was too big, division by Inf has given NA - convert to 1s
prob_lst <- rapply(prob_lst, function(x){ifelse(is.na(x), 1, x)},
                    how = "replace")


# create output directory
#dirname <- paste0("output/plots/", ID, "_", gsub(":", "-", gsub(" ", "T", Sys.time())))
#dir.create(dirname, recursive = TRUE)


# plotting loop
for(i in 1:length(prob_lst)){

# create colour palette
col.pal <- hcl.colors(length(prob_lst[[i]]), palette = "viridis")

# open plotting connection
jpeg(paste0(dirname, "/", names(prob_lst)[i], ".jpg"),
     width = 85, height = 85, units = 'mm', res = 300)

# base plot to provide axis labels, actual lines are omitted
matplot(prob_lst[[1]][[1]], type="n",
        ylab="Probability of sharing pest species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F)

# add one line per model
for(j in 1:length(prob_lst[[i]])){
lines(prob_lst[[i]][[j]], col=col.pal[j])
}

# legend
for(j in 1:length(prob_lst[[i]])){
if(j==1){
leg <- paste0("Model ", j)
}else{
leg <- c(leg, paste0("Model ", j))
}}

legend("topright",
       legend = leg,
       col = col.pal,
       lty = c(1, 1),
       lwd = c(3, 3),
       cex = 0.3,
       bty = "n")

# x axis (below)
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis = "black", las = 1, cex.axis = 0.5)
# y axis (left)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis = "black", las = 1, cex.axis = 0.5)


# close connection
dev.off()
}

