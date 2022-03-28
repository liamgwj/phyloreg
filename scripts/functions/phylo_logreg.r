# LJ 2022-02-26
# perform phylogenetic logistic regression

# requires objects 'phy', 'db', and 'runID'

library(ape)

# prep data for model ---------------------------------------------------------

# get all pairwise phylogenetic distances
pd <- cophenetic.phylo(phy)

# log transform
logpd <- log10(pd + 1)

# get vector of all pest names
pests <- rownames(db)

# get vector of all host names
hosts <- colnames(db)

# prep output data.frame
coef <- data.frame()

# logistic regression model proper --------------------------------------------

# loop through all potential pests --------------------------------
for(i in seq_along(pests)){

# select focal pest for this iteration
fpest <- pests[i]

# get all hosts of focal pest
host_subset <- hosts[which(db[fpest, ]>0)]

# loop through all known hosts of focal pest -----------------------------
for(j in seq_along(host_subset)){

# select focal host for this sub-iteration
fhost <- host_subset[j]

# subset data for model run and make sure observation order matches
compatibility <- as.numeric(as.vector(db[fpest, ][sort(names(db[fpest, ]))]))
logpd_subset <- logpd[, fhost][sort(names(logpd[, fhost]))]

# fit logistic regression model
## should specify na.action
lrm_out <- stats::glm(compatibility ~ logpd_subset,
                      family = stats::binomial(link = "logit"))

# organize model coefficients
stats_lrm_out <- stats::coef(summary(lrm_out))

conf_int <- stats::confint.default(lrm_out)

run_coef <- data.frame(fpest,
                       fhost,
                       lrm_out$coefficients[1],
                       stats_lrm_out[1, 2],
                       stats_lrm_out[1, 3],
                       stats_lrm_out[1, 4],
                       conf_int[1, 1],
                       conf_int[1, 2],
                       lrm_out$coefficients[2],
                       stats_lrm_out[2, 2],
                       stats_lrm_out[2, 3],
                       stats_lrm_out[2, 4],
                       conf_int[2, 1],
                       conf_int[2, 2])

colnames(run_coef) <- c("Pest", "FocalHost", "intercept", "Std. Error",
                        "z value", "Pr(>|z|)", "2.5 %", "97.5 %", "slope",
                        "Std. Error", "z value", "Pr(>|z|)", "2.5 %", "97.5 %")

rownames(run_coef) <- NULL

# append to output data.frame
coef <- rbind(coef, run_coef)

} # end of host loop (j)
} # end of pest loop (i) ---------------------------------------------

# write out -------------------------------------------------------------------


if(!dir.exists("output")){
dir.create("output")
}

# all coefficients
write.csv(coef,
          paste0("output/", runID, subID, "_coef.csv"),
          row.names = FALSE)

# average coefficients for each pest
coef_av <- aggregate(coef[,-which(colnames(coef)%in%c("Pest", "FocalHost"))],
                     list(Pest = coef$Pest),
                     mean)

write.csv(coef_av,
          paste0("output/", runID, subID, "_coef-av.csv"),
          row.names = FALSE)

