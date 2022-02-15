
setwd("/home/liam/Documents/MSc/analysis/sensitivity_analysis/simulated")

# LJ 2022-02-08 simulate a simple phylogeny and vary placement of "host" tips

# load packages
library(ape)
library(geiger)
library(TreeSim)
library(geotax)
library(dplyr)

# set date/time ID
now <- gsub(" ", "T", Sys.time())


# simulate phylogenies
phy <- sim.bd.taxa(n = 300,
                   numbsim = 3,
                   lambda = 0.3,
                   mu = 0.25,
                   complete = FALSE,
                   stochsampling = TRUE,
                   frac = 0.75
                   )


# create transition matrix
q <- list(rbind(c(-.01, .01), # probability of state change 1 -> 2 
                c(.05, -.05) # probablility of state change 2 -> 1 (reversion)
                ))

# simulate trait evolution on phylogenies
char <- lapply(phy, function(x){sim.char(phy = x,
                                         par = q,
                                         nsim = 1,
                                         model = "discrete",
                                         root = 1
                                         )})


# identify compatible tips
hosts <- lapply(char, function(x){rownames(x)[which(x==2)]})


nonhosts <- lapply(char, function(x){rownames(x)[which(x==1)]})


for(k in 1:3){

if(k==1){frac <- 1}
if(k==2){frac <- 0.5}
if(k==3){frac <- 0.25}

tokeep <- list()

todrop <- list()

for(l in 1:length(nonhosts)){
tokeep[[l]] <- sample(nonhosts[[l]], length(nonhosts[[l]])*frac)

todrop[[l]] <- setdiff(nonhosts[[l]], tokeep[[l]])

phy[[l]] <- drop.tip(phy[[l]], todrop[[l]])
}



# get all pairwise phylogenetic distances
pd <- lapply(phy, cophenetic.phylo)

# log transform pd values
logpd <- lapply(pd, function(x){log10(x + 1)})





for(j in seq_along(hosts)){

out <- data.frame()

out_data_all <- data.frame(pd=seq(0, 650, 1))

for(i in seq_along(hosts[[j]])){

# choose focal host
focal <- hosts[[j]][i]

# recode and rearrange character info
char.bin <- char[[j]]

if(frac<1){
char.bin <- char.bin[-which(!rownames(char.bin)%in%phy[[j]]$tip.label)]
}

char.bin[which(char.bin==1)] <- 0
char.bin[which(char.bin==2)] <- 1
char.bin <- matrix(char.bin, nrow=1)

if(frac<1){
colnames(char.bin) <- rownames(char[[j]])[-which(!rownames(char[[j]])%in%phy[[j]]$tip.label)]
}else{
colnames(char.bin) <- rownames(char[[j]])
}


run.pd <- logpd[[j]][focal, ]

run.dat <- rbind(char.bin, run.pd)

# fit logistic regression model
log_out <- stats::glm(run.dat[1,] ~ run.dat[2,],
                          family = stats::binomial(link = "logit"))


    stats_log_out <- stats::coef(summary(log_out))
    
    conf_int <- stats::confint.default(log_out)
    
    lrcoeffs <- data.frame(log_out$coefficients[1],
                           stats_log_out[1, 2],
                           stats_log_out[1, 3],
                           stats_log_out[1, 4],
                           conf_int[1, 1],
                           conf_int[1, 2],
                           log_out$coefficients[2],
                           stats_log_out[2, 2],
                           stats_log_out[2, 3],
                           stats_log_out[2, 4],
                           conf_int[2, 1],
                           conf_int[2, 2])
    
    colnames(lrcoeffs) <- c("intercept", "Std. Error", "z value", "Pr(>|z|)",
                            "2.5 %", "97.5 %", "slope", "Std. Error",
                            "z value", "Pr(>|z|)", "2.5 %", "97.5 %")
    
    rownames(lrcoeffs) <- NULL
    


coef <- lrcoeffs


## which function is this from? seems to have been multiplying SD (coef[2])
# predict suitability
logit <- coef[[1]] + coef[[7]] * run.pd

prob_logit <- (exp(logit)/(1 + exp(logit)))

tmp <- data.frame(t(prob_logit))

tmp$phy <- j

tmp$focal <- focal

out <- rbind(out, tmp)


# write out data for combined plotting later
# predicted suitability for each pest at PDs 0 - 650 MY
out_data <- apply(coef[ , c(1, 7)],
                  1,
                  function(x) prob_logit(x, log10(seq(0, 650, 1) + 1))
                  ) %>%
            as.data.frame()

colnames(out_data) <- hosts[[j]][i]

out_data_all <- cbind(out_data_all, out_data)


}

ntip <- Ntip(phy[[j]])

out2 <- out[,c(ntip+1,ntip+2,1:ntip)]


write.csv(out2, paste0("output/psuit_phy", j, "_", frac, "_", now, ".csv"), row.names=FALSE)


write.csv(out_data_all, paste0("output/curvedata", k,j,".csv"), row.names = FALSE)

}

write.tree(phy, paste0("output/multiphy_", frac, "_", now, ".nwk"))
}


