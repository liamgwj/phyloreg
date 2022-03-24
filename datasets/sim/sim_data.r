# LJ 2022-02-08 simulate a simple phylogeny and "host" tips

setwd("/home/liam/Documents/MSc/analysis/phyloreg/datasets/sim/")

# load packages
library(ape)
library(geiger)
library(TreeSim)

# set date/time ID
now <- gsub("[:-]", "", gsub(" ", "T", Sys.time()))


# simulate phylogeny/ies
phy <- sim.bd.taxa(n = 20,
                   numbsim = 1,
                   lambda = 0.3,
                   mu = 0.25,
                   complete = FALSE,
                   stochsampling = TRUE,
                   frac = 0.75
                   )


# create transition matrix
q <- list(rbind(c(-.1, .1), # probability of state change 1 -> 2 
                c(.05, -.05) # probablility of state change 2 -> 1 (reversion)
                ))

# simulate trait evolution on phylogeny/ies
char <- lapply(phy, function(x){sim.char(phy = x,
                                         par = q,
                                         nsim = 1,
                                         model = "discrete",
                                         root = 1
                                         )})


# write out

if(!dir.exists("simulated_data")){
dir.create("simulated_data")
}

write.tree(phy,
           paste0("simulated_data/sim-", now, "_phy.nwk"))

write.csv(ifelse(t(data.frame(char[[1]]))==2, 1, 0),
          paste0("simulated_data/sim-", now, "_db.csv"))

