# LJ 2022-03-09 plot model predictions presentably on phylogeny

# add line type - dashed lines for tips not used in model fitting

# setup -----------------------------------------------------------------------

library(ape)
library(viridis)

setwd("/home/liam/Documents/MSc/analysis/phyloreg/LJ_Parker")

# read in phylogeny
# modified QJ/Z (missing 2 genera/tips)
phy <- read.tree("output/QJ_Parker_2022-03-09.nwk")

# Parker original, with my additions (many polytomies)
#phy <- read.tree("output/cleaned_Parker_phylogeny_2022-02-28.nwk")


# read in model coefficients
coef_1 <- read.csv("output/model_coefficients/parker-QJ-all_coef-av_2022-03-09T12-31-19.csv")
coef_2 <- read.csv("output/model_coefficients/parker-QJ-halfhosts_coef-av_2022-03-10T15-38-24.csv")


# read in host/pest db
db <- read.csv("output/cleaned_Parker_database_2022-02-28.csv",
               row.names = 1)


# generate tip-level predictions for plotting ---------------------------------

# format list of coefficients
coef_lst <- vector("list", length = nrow(coef_2))

names(coef_lst) <- coef_2$Pest

for(i in 1:nrow(coef_2)){
coef_lst[[i]] <- rbind(coef_1[which(coef_1$Pest == coef_2$Pest[i]) ,],
                       coef_2[i,])
}

# choose a pest
i=1
fpest <- coef_lst[[i]]$Pest[1]

fp_df <- data.frame(t(db[fpest,]))
colnames(fp_df) <- "compatible"
fp_df$avPD <- NA

fp_df <- fp_df[which(rownames(fp_df)%in%phy$tip.label),]

hosts <- rownames(fp_df)[fp_df$compatible==1]

for(i in seq_along(hosts)){
fp_df[,ncol(fp_df)+1] <- cophenetic(phy)[,hosts[i]]
colnames(fp_df)[ncol(fp_df)] <- paste0("PDfh", i)
}

fp_df$avPD <- rowMeans(fp_df[,(ncol(fp_df)-length(hosts)+1):ncol(fp_df)])

fp_df$minPD <- apply(fp_df[,(ncol(fp_df)-length(hosts)+1):ncol(fp_df)],1,min)



fp_df$logit1min <- coef_lst[[i]]$intercept[1] + coef_lst[[i]]$slope[1] * log10(fp_df$minPD + 1)

fp_df$pred1min <- exp(fp_df$logit1min)/(1 + exp(fp_df$logit1min))

# if predicted value was too big, division by Inf has given NA - convert to 1s
fp_df$pred1min <- ifelse(is.na(fp_df$pred1min), 1, fp_df$pred1min)

i=1

fp_df$logit2min <- coef_lst[[i]]$intercept[2] + coef_lst[[i]]$slope[2] * log10(fp_df$minPD + 1)

fp_df$pred2min <- exp(fp_df$logit2min)/(1 + exp(fp_df$logit2min))

# if predicted value was too big, division by Inf has given NA - convert to 1s
fp_df$pred2min <- ifelse(is.na(fp_df$pred2min), 1, fp_df$pred2min)



# plotting palette
#cols <- rev(heat.colors(101))
cols <- viridis(101, option="inferno")

# give terminal edges a colour dependent on their predicted suitability value

edgecols1 <- rep("black", Nedge(phy))

for(i in seq_along(rownames(fp_df))){
edgecols1[which.edge(phy, rownames(fp_df)[i])] <- cols[(round(fp_df$pred1min[i], 2)*100)+1]
}

for(i in seq_along(hosts)){
edgecols1[which.edge(phy, hosts[i])] <- "blue"
}


edgecols2 <- rep("black", Nedge(phy))

for(i in seq_along(rownames(fp_df))){
edgecols2[which.edge(phy, rownames(fp_df)[i])] <- cols[(round(fp_df$pred2min[i], 2)*100)+1]
}

for(i in seq_along(hosts)){
edgecols2[which.edge(phy, hosts[i])] <- "blue"
}



# plot ------------------------------------------------------------------------

depth <- max(node.depth.edgelength(phy))

par(mfrow=c(1,2))
plot(phy,
     x.lim = depth+(depth/10),
     label.offset = 5,
     no.margin = TRUE,
     edge.color = edgecols1)
plot(phy,
     direction = "leftwards",
     show.tip.label = FALSE,
     edge.color = edgecols2)




