# LJ 2022-02-11 plot output from simulated logreg test

setwd("/home/liam/Documents/MSc/analysis/sensitivity_analysis/simulated/output")

library(ape)

# read in sim output
psuit1 <- read.csv("psuit_phy1_1_2022-02-14T12:23:31.csv")
psuit5 <- read.csv("psuit_phy1_0.5_2022-02-14T12:23:31.csv")
psuit25 <- read.csv("psuit_phy1_0.25_2022-02-14T12:23:31.csv")

mphy <- read.tree("multiphy_1_2022-02-14T12:23:31.nwk")

# plotting palette
cols <- heat.colors(11)

for(i in 1:3){

if(i==1){ d <- psuit1}
if(i==2){ d <- psuit5}
if(i==3){ d <- psuit25}

hosts <- unique(d$focal)

av <- colMeans(d[3:ncol(d)])

avRound <- round(av, 1)

#hist(avRound)
#bootstrapping vs. "dynamic iteration" might make a difference here?

edge_colours <- rep("black", Nedge(mphy[1]))


# give terminal edges a colour dependent on their predicted suitability value
for(j in 1:Ntip(mphy[1])){

if(mphy[[1]]$tip.label[j]%in%names(avRound)){

        edge_colours[which.edge(mphy[[1]], mphy[[1]]$tip.label[j])] <-

                cols[11-(avRound[which(names(avRound) == mphy[[1]]$tip.label[j])] *10)]

}}

tipcol <- data.frame(col=rep("black", Ntip(mphy[1])))
tipcol$tip <- mphy[[1]]$tip.label
tipcol$col[which(tipcol$tip%in%hosts)] <- "red"

# plot
png(paste0("outphy_", i, ".png"),  width = 600, height = 600, units = "px")
plot(mphy[[1]],
     edge.color = edge_colours,
     tip.color = tipcol$col[match(tipcol$tip,
                                      mphy[[1]]$tip.label)]
)
dev.off()

png(paste0("hist_", i, ".png"),  width = 600, height = 600, units = "px")
hist(av)
dev.off()
}



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
ps1_3 <- read.csv("curvedata11.csv")
ps1_10 <- read.csv("curvedata21.csv")
ps1_22 <- read.csv("curvedata31.csv")


ps1_3$av <- rowMeans(ps1_3[,2:ncol(ps1_3)])
ps1_10$av <- rowMeans(ps1_10[,2:ncol(ps1_10)])
ps1_22$av <- rowMeans(ps1_22[,2:ncol(ps1_22)])


# create plot and write to file
jpeg(paste0("averaged.jpg"), width = 85, height = 85, units = 'mm', res = 300)

# base plot to provide axis labels, actual lines are omitted
matplot(ps1_3[,2:ncol(ps1_3)], type="n",
        ylab="Probability of sharing pest species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F
        )

# one line per model run
lines(ps1_3[,24], col = bluish_green, lwd = 3, lty = 1)
lines(ps1_10[,24], col = sky_blue, lwd = 3, lty = 1)
lines(ps1_22[,24], col = orange, lwd = 3, lty = 1)

# legend
legend(300, 1,
       legend = c("whole tree",
                  "50% nonhost pruned",
                  "75% nonhost pruned"),
       lty = c(1, 1, 1),
       col = c(bluish_green, sky_blue, orange),
       lwd = c(3, 3, 3),
       cex = 0.3,
       bty = "n")


abline(v=max(node.depth.edgelength(mphy[[1]]))*2)


# x axis (below)
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis = "black", las = 1, cex.axis = 0.5)
# y axis (left)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis = "black", las = 1, cex.axis = 0.5)

dev.off()



