
# this workflow is a little silly - using averaged model coefficients to make
# per-focal host predictions, then averaging those predictions again...

# Q: better to bootstrap with large n or have exactly one run per focal host?

# setup -----------------------------------------------------------------------

# devtools::install_github("alrobles/geotax")
require(geotax)
library(magrittr)
library(ape)

# host phylogenetic tree built with phylomatic
phy_host <- read.tree("indata/host.tree")

# full beetle/host tree compatibility database
db_pest <- read.csv("indata/full_taxa_beetle_plant.csv", header = T)

# expand database to matrix and remove columns (host plants) not in phy
mx_pest <- incidence(db_pest[ , c(7,13)])[ , phy_host$tip.label]

# subset to only beetle species having > 2 hosts
mx_pest <- mx_pest[rowSums(mx_pest) > 2, ]


# model loop ------------------------------------------------------------------

for(p in seq_along(rownames(mx_pest))){

for(k in 1:3){

# set number of non-focal pests to include
if(k==1){n <- 3}
if(k==2){n <- 10}
if(k==3){n <- 19}

# check for output directories and create if necessary
if(!dir.exists("output")){dir.create("output")}
if(!dir.exists(paste0("output/", rownames(mx_pest)[p]))){
        dir.create(paste0("output/", rownames(mx_pest)[p]))}
if(!dir.exists(paste0("output/", rownames(mx_pest)[p], "/focalPlus_", n, "_pests"))){
        dir.create(paste0("output/", rownames(mx_pest)[p], "/focalPlus_", n, "_pests"))}


# subset to one focal pest plus n random others
ind_pest_focal <- which(rownames(mx_pest) == rownames(mx_pest)[p])

ind_pest_all <- 1:length(rownames(mx_pest))

mx_pest_subset <- mx_pest[c(ind_pest_focal,
                            sample(ind_pest_all[-ind_pest_focal],
                                   n, replace = FALSE)
                            ), ]

# remove host plants without included pests
mx_pest_subset <- mx_pest_subset[ , which(colSums(mx_pest_subset) > 0)]


# prune phy to only included host species
phy_host_subset <- keep.tip(phy_host, colnames(mx_pest_subset))

# get all pairwise phylogenetic distances
pd_host_subset <- phy_host_subset %>% cophenetic.phylo()

# for plotting later, also get ALL PD values
pd_host <- phy_host %>% cophenetic.phylo()

# log transform pd values
logpd_host_subset <- log10(pd_host_subset + 1)


# calculate logistic regression coefficients for each species
coef_pest_subset <- sapply(1:nrow(mx_pest_subset),
                   function(x) log_reg_boostrap(mx_pest_subset[x, , drop = F],
                                                logpd_host_subset,
                                                1000)
                           ) %>%
                    t()


# write out data for combined plotting later
# predicted suitability for each pest at PDs 0 - 650 MY
out_data <- apply(coef_pest_subset[ , c(1, 7)],
                  1,
                  function(x) prob_logit(x, log10(seq(0, 650, 1) + 1))
                  ) %>%
            as.data.frame()

colnames(out_data) <- rownames(mx_pest_subset)

write.csv(out_data, paste0("output/", rownames(mx_pest)[p], "/focalPlus_", n, "_pests/curve_data.csv"), row.names = FALSE)


## phylogeny plots

# for focal pest, loop through known hosts and plot predicted suitability for
# each tip on the phylogeny

# all known hosts of focal pest
hosts_pest_focal <- names(mx_pest[ind_pest_focal,
                                  which(mx_pest[ind_pest_focal, ] == 1)])

# df to hold averaged predictions
pred_av <- data.frame(host_gen = names(pd_host[1, ]))

for(i in seq_along(hosts_pest_focal)){

# isolate focal host for run
focal_host <- hosts_pest_focal[i]

# get pairwise PD values from focal host to all other hosts
pred_i <- data.frame(pd_host[focal_host, ])

# get modelled probability of suitability for all hosts
pred_i$pred <- prob_logit(coef_pest_subset[1, c(1, 7)], log10(pred_i[,1] + 1))

# round to 1 decimal place for plotting
pred_i$pred_round <- round(pred_i$pred, 1)

# copy values for later averaging
pred_av[ , (i+1)] <- pred_i$pred

# create palette of 11 colours
cols <- heat.colors(11)

# create edge colour object, default colour is black
# note that this includes colours for internal branches
edge_colours <- rep("black", Nedge(phy_host))

# give terminal edges a colour dependent on their predicted suitability value
for(j in 1:Ntip(phy_host)){
        
        edge_colours[which.edge(phy_host,
                                phy_host$tip.label[j])] <- 
                
                cols[11-(
                        
                        pred_i[which(rownames(pred_i) == 
                                             phy_host$tip.label[j]),3]
                        
                        *10)]

        }

# colour edge leading to focal host blue
edge_colours[which.edge(phy_host, focal_host)] <- "blue"



# colour tip names of known hosts red
pred_i$namecol <- "black"
pred_i$namecol[which(rownames(pred_i) %in% hosts_pest_focal)] <- "red"

# plot
png(paste0("output/", rownames(mx_pest)[p], "/focalPlus_", n, "_pests/focalHost_", i, ".png"))
plot(phy_host,
     edge.color = edge_colours,
     tip.color = pred_i$namecol[match(rownames(pred_i),
                                      phy_host$tip.label)]
        )
dev.off()


} # end of phylogeny plotting loop


## average predictions plots

pred_av$av <- rowMeans(pred_av[ , 2:ncol(pred_av)])

# round to 1 decimal place for plotting
pred_av$pred_round <- round(pred_av$av, 1)

# create edge colour object, default colour is black
# note that this includes colours for internal branches
edge_colours <- rep("black", Nedge(phy_host))


# give terminal edges a colour dependent on their predicted suitability value
for(j in 1:Ntip(phy_host)){
        
        edge_colours[which.edge(phy_host, phy_host$tip.label[j])] <- 
                
                cols[11-(pred_av[which(pred_av[, 1] == phy_host$tip.label[j]), "pred_round"] *10)]
        
}


# plot
png(paste0("output/", rownames(mx_pest)[p], "/focalPlus_", n, "_pests/averaged.png"))
plot(phy_host,
     edge.color = edge_colours,
     tip.color = pred_i$namecol[match(rownames(pred_i),
                                      phy_host$tip.label)]
)

legend(10, 10,
       legend = c("1",
                  "0.5",
                  "0.1"),
       lty = c(1, 1, 1),
       col = c(cols[1], cols[6], cols[10]),
       lwd = c(3, 3, 3),
       cex = 0.5,
       bty = "n")

dev.off()


} # end of model loop


# combined model output plot --------------------------------------------------

# set colour palette
orange <- rgb(red = 230, green = 159, blue = 0, max = 255)
sky_blue <- rgb(86, 180, 233, max = 255)
bluish_green <- rgb(0, 158, 115, max = 255)
vermillion <- rgb(213, 94, 0, max = 255)
reddish_purple <- rgb(204, 121, 167, max = 255)
yellow <- rgb(240, 228, 66, max = 255)
white <- rgb(255, 255, 255, max = 255)

# read in model loop output
## names need updating
ps1_3 <- read.csv(paste0("output/", rownames(mx_pest)[p], "/focalPlus_3_pests/curve_data.csv"))
ps1_10 <- read.csv(paste0("output/", rownames(mx_pest)[p], "/focalPlus_10_pests/curve_data.csv"))
ps1_22 <- read.csv(paste0("output/", rownames(mx_pest)[p], "/focalPlus_19_pests/curve_data.csv"))

# create plot and write to file
jpeg(paste0("output/", rownames(mx_pest)[p], "/curves.jpg"), width = 85, height = 85, units = 'mm', res = 300)

# base plot to provide axis labels, actual lines are omitted
matplot(ps1_3, type="n",
        ylab="Probability of sharing beetle species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F
        )

# one line per model run
lines(ps1_3[,1], col = bluish_green, lwd = 3, lty = 1)
lines(ps1_10[,1], col = sky_blue, lwd = 3, lty = 1)
lines(ps1_22[,1], col = orange, lwd = 3, lty = 1)

# legend
legend(300, 1,
       legend = c("Focal + 3",
                  "Focal + 10",
                  "Focal + 19"),
       lty = c(1, 1, 1),
       col = c(bluish_green, sky_blue, orange),
       lwd = c(3, 3, 3),
       cex = 0.3,
       bty = "n")

# x axis (below)
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis = "black", las = 1, cex.axis = 0.5)
# y axis (left)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis = "black", las = 1, cex.axis = 0.5)

dev.off()

}
