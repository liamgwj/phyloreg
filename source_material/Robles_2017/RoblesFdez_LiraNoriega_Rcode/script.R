#Support material
#path the the material directory
setwd("C:/Users/Andres/Desktop/material")

devtools::install_github("alrobles/geotax")
require(geotax)


install.packages("bossMaps")
install.packages("magrittr")
install.packages("broom")
install.packages("biglm")
install.packages("ape")
install.packages("rasterVis")
install.packages("raster")

pack <- c("magrittr", "broom", "biglm", "ape", "bossMaps", "rasterVis", "raster")
lapply(pack, require, character.only = TRUE)
  
####host phylogenetical tree build with phylomatic
tree <- read.tree("host.tree")

#distance matrix
d <- tree %>% cophenetic.phylo()
d <- log10(d+1)

#not run
#d <- log10(phy_dist + 1) # using data object

#############################################
# create incidence matrix from database     #
# and filter with distance matrix           #
#############################################

db <- read.csv("full_taxa_beetle_plant.csv", header = T)
i_matrix <- incidence(db[, c(7, 13)])[ ,colnames(d)]


# filtering incidence matrix with incidence > 2
i <- i_matrix[rowSums(i_matrix)>2, ]

# distance matrix

#logistic regression coefficients for all the set
coef <- log_reg_boostrap(i, d, 1000)

#logistic regression coefficients for each case
coef_all <- sapply(1:nrow(i), function(x) log_reg_boostrap(i[x, ,drop=F],
                                                           d, 1000) ) %>% t
#filter matrix for each case taking only source infected host to all target host
d_all <- lapply(1:nrow(i), function(x){ d[which(i[x, ]  %in% 1), ] })
names(d_all) <-  rownames(i)

#probability matrix for all set
p <- prob_logit(coef[c(1,7)], d)

#probability matrix for each case
p_all <- lapply(1:nrow(coef_all),
                function(x) prob_logit(coef_all[x, c(1,7)], d_all[[x]]) )
names(p_all) <- rownames(i)

m1_db <- read.csv("m1_pts.csv") #points database retrived from gbif by species
m2_db <- read.csv("m2_pts.csv") #points database retrived from gbif by genera

M1 <-  read.csv("M1.csv")   #presense-absence matrix from m1_db at 1 degree resolution
row.names(M1) <-  M1[ ,1]
M1 <- M1[ ,-1]

M2 <-  read.csv("M2.csv")   ##presense-absence matrix from m1_db 1 degree resolution
row.names(M2) <-  M2[ ,1]
M2 <- M2[ ,-1]

#M1 and M2 can be build with geotax::PAM() with geotax::world shape

#match data from genera in points databases and probability matrix

n1 <- colnames(p)[ unique(m1_db[ ,1]) %>% as.matrix %>%
                     match(., colnames(p)) %>% na.omit ]
n2 <- colnames(p)[ unique(m2_db[ ,1]) %>% as.matrix %>%
                     match(., colnames(p)) %>% na.omit ]

M1 <- M1[n1, ] %>% as.matrix  #presence-absence matrix M1 filter with match
M2 <- M2[n2, ] %>% as.matrix  #presence-absence matrix M2 filter with match

p_all_m1 <- lapply(p_all, function(x) x[, n1] ) #probability matrix filter
p_all_m2 <- lapply(p_all, function(x) x[, n2] ) #probability matrix with filter

p_m1 <- p[n1, n1]
p_m2 <- p[n2, n2]

#calculate geographical phylogenetical information

#full set case
G1 <- p_m1 %*% M1
G2 <- p_m2 %*% M2

#all single case
G1_all <- lapply(p_all_m1, function(x) x %*% M1) #G matrix
G2_all <- lapply(p_all_m2, function(x) x %*% M2)

g1 <- richness_PAM(world, G1, 1 ) #g raster
g1_log <- 1 - g1 %>% normalize() %>% cumulative() %>% logistic()

g2 <- richness_PAM(world, G2, 1 )
g2_log <- 1 - g2 %>% normalize() %>% cumulative() %>% logistic()

g1_all <-  lapply(G1_all, function(x) richness_PAM(world, x, 1 ) )
g2_all <-  lapply(G2_all, function(x) richness_PAM(world, x, 1 ) )

g1_log_all <- lapply(g1_all, function(x){
  1 - x %>% normalize() %>% cumulative() %>% logistic() } )
g2_log_all <- lapply(g2_all, function(x){
  1 - x %>% normalize() %>% cumulative() %>% logistic() } )

richness_M1 <- richness_PAM(world, M1, 1 )
richness_M2 <- richness_PAM(world, M2, 1 )

cor_r_g1_M1 <- lapply(g1_log_all, function(i){
  x <- i@data@values
  y <- richness_M1@data@values
  cor(x,y)
}) %>% do.call(rbind, .)

cor_r_g2_M2 <- lapply(g2_log_all, function(i){
  x <- i@data@values
  y <- richness_M2@data@values
  cor(x,y) }) %>% do.call(rbind, .)

r_vs_h <- data.frame(rowSums(i), exp(-rowSums(i)),
                     cor_r_g1_M1, cor_r_g2_M2)

colnames(r_vs_h) <-  c("host_range", "h", "r_M1", "r_M2" )

model_M1 <-  lm(r_M1 ~ h , r_vs_h)
model_M2 <-  lm(r_M2 ~ h , r_vs_h)


###
# graphics
###

###color pallete

orange <- rgb(red= 230,green= 159,blue= 0, max = 255)
sky_blue <- rgb(	86,180,233, max = 255)
bluish_green	<- rgb(0,158,115, max = 255)
vermillion <- rgb(213,94,0, max = 255)
reddish_purple <- rgb(204,121,167, max = 255)
yellow <-  rgb (240,228,66, max = 255)
white <-  rgb(255,255,255, max = 255 )

####
# probability phylogenetic distance
PD    <- seq(0, 650, 1) #x axis values across range of phylogenetic distances
ps1   <- apply(coef_all [ ,c(1,7)], 1, function(x) prob_logit(x, log10(PD+1) ) ) %>% as.data.frame()
colnames(ps1) <-  rownames(i)
ps_all   <- prob_logit(coef[c(1,7)], log10( PD + 1) )


jpeg("1.jpg", width = 85, height = 85, units = 'mm', res = 300)
matplot(ps1, type="l",
        ylab="Probability of sharing beetle species",
        xlab="Phylogenetic distance from source to target host (My)",
        cex.lab = 0.5,
        cex=0.5, axes = F,
        col="gray",lwd=0.1,lty=1 )

lines(ps1$`Xyleborus xylographus`, col = yellow, lwd=3, lty = 1)
lines(ps1$`Xyleborus glabratus`, col = bluish_green, lwd=3, lty = 1 )
lines(ps1$`Xylosandrus crassiusculus`, col = sky_blue, lwd=3, lty = 1 )
lines(ps_all, col = vermillion, lwd=2, lty =2 )

legend(300, 1,
       legend = c("Each beetle","", "Xyleborus xylographus", "(Narrow range, phylogenetically constrained)",
                  "Xyleborus glabratus", "(Wide range, phylogenetically constrained)",
                  "Xylosandrus crassiusculus", "(Wide range with phylogenetical signal)",
                  "All beetles"),
       lty= c(1,0,1,0,1,0,1,0,2), col = c("grey", "", yellow,"", bluish_green,
                                          "",sky_blue,"", vermillion ),
       lwd= c(2,0,3,0,3,0,3,0,2), cex = 0.3, bty="n")
axis(1, seq(0, 650, 50), seq(0, 650, 50),
     col.axis= "black", las=1, cex.axis=0.5)
axis(2, seq(0, 1, 0.1), seq(0, 1, 0.1),
     col.axis= "black", las=1, cex.axis=0.5)
dev.off()

################
# raster plots #
################

all <- mask(stack(g1_log, g2_log), world)
names(all) <- c("M1", "M2")

cases_M1 <- stack(g1_log_all$`Xyleborus xylographus`,
                  g1_log_all$`Xyleborus glabratus`,
                  g1_log_all$`Xylosandrus crassiusculus`) %>% mask(., world)
cases_M2 <- stack(g2_log_all$`Xyleborus xylographus`,
                  g2_log_all$`Xyleborus glabratus`,
                  g2_log_all$`Xylosandrus crassiusculus`) %>% mask(.,world)
cases <- stack(cases_M1, cases_M2)

#PLOT THEME
myTheme <- BTCTheme()
myTheme$panel.background$col = 'gray'
rainbTheme5 <- rasterTheme(region = rev(rainbow(n = 5)))

#PLOT

jpeg("2.jpg", width = 180, height = 85, units = 'mm', res = 300)
levelplot(all, colorkey=list(space="bottom"),
          par.settings= myTheme,
          names.attr = rep(" ", 2 ) )
dev.off()

jpeg("2_bin.jpg", width = 180, height = 85, units = 'mm', res = 300)
levelplot(all > 0.5, colorkey=F,
          par.settings= myTheme,
          names.attr = rep(" ", 2 ) )
dev.off()


jpeg("3.jpg", width = 180, height = 180, units = 'mm', res = 300)
levelplot(cases, colorkey=list(space="bottom"),
          par.settings = myTheme,
          names.attr = rep(" ", nlayers(cases)) )
dev.off()

jpeg("3_bin.jpg", width = 180, height = 180, units = 'mm', res = 300)
levelplot(cases>0.5, colorkey=F,
          par.settings= myTheme,
          names.attr = letters[1:6] )
dev.off()

###
# plot regression with richness
##

eq_M1 <- substitute(italic(r) == a + b %.% italic(h)*","~~~italic(r)^2~"="~r2*","~~~italic("p-value")~"="~p,
                    list(a = format(coef(model_M1)[1], digits = 3),
                         b = format(coef(model_M1)[2], digits = 3),
                         r2 = format(summary(model_M1)$r.squared, digits = 3),
                         p = format(glance(model_M1)$p.value)))
f_M1 <- as.expression(eq_M1)
eq_M2 <- substitute(italic(r) == a + b %.% italic(h)*","~~~italic(r)^2~"="~r2*","~~~italic("p-value")~"="~p,
                    list(a = format(coef(model_M2)[1], digits = 3),
                         b = format(coef(model_M2)[2], digits = 3),
                         r2 = format(summary(model_M2)$r.squared, digits = 3),
                         p = format(glance(model_M2)$p.value)))
f_M2 <- as.expression(eq_M2)


jpeg("4.jpg", width = 85, height = 85, units = 'mm', res = 300)
plot(NULL, xlim=c(0, max(r_vs_h$h)), ylim= c(0.5,1),
     ylab="Richness correlation", xlab= expression('h = e'^-'H'),
     cex.axis = 0.5, cex.lab = 0.5)
abline(model_M1, col = vermillion, lwd= 2)
abline(model_M2, col = bluish_green, lwd= 2)
points(r_M1 ~ h, data = r_vs_h,  pch=21,  bg = vermillion, cex = 0.5)
points(r_M2 ~ h, data = r_vs_h,  pch=21,  bg = bluish_green, cex = 0.5)

legend("topright",legend=c("M1", "M2"),
       col=c(vermillion, bluish_green),
       pch= c(16, 16), bty="n",
       ncol=1, cex=0.5, pt.cex=0.5)

legend("bottomleft", legend=c(f_M1, f_M2),
       col = c(vermillion, bluish_green),
       pch = c(16, 16), bty = "n",
       ncol= 1, cex = 0.4, pt.cex = 0.4,
       text.font= 2 )


write.csv(round(coef_all[, c(1,5,6,7, 11,12)], 4), "tabla_1.csv")
rownames(coef_all)  <-  rownames(i)

