# LJ 2022-02-08 simulate a simple phylogeny and "host" tips

setwd("/home/liam/Documents/MSc/analysis/phyloreg/r")

# load packages
library(ape)
library(geiger)
library(TreeSim)

# set date/time ID
now <- gsub("[:-]", "", gsub(" ", "T", Sys.time()))

# set number of pests to simulate on single phylogeny
npest <- 3

# simulate phylogeny/ies
phy <- sim.bd.taxa(n = 5,
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

# simulate trait evolution on phylogeny
for(i in 1:npest){

char <- lapply(phy, function(x){sim.char(phy = x,
                                         par = q,
                                         nsim = 1,
                                         model = "discrete",
                                         root = 1
                                         )})

row_i <- ifelse(t(data.frame(char[[1]]))==2, 1, 0)

if(i=1){
out <- row_i}else{
out <- rbind(out, row_i)}
}

# write out
if(!dir.exists("indata/simulated")){
dir.create("indata/simulated")
}

write.tree(phy,
           paste0("indata/simulated/sim-", now, "_phy.nwk"))

write.csv(out,
          paste0("indata/simulated/sim-", now, "_db.csv"))

