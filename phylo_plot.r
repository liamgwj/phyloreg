# LJ 2022-03-09 plot model predictions presentably on phylogeny

# something's still up - all sisters in a polytomy should be same p/colour

# setup -----------------------------------------------------------------------

library(ape)
library(viridis)


# generate tip-level predictions for plotting ---------------------------------

# choose a pest
j=2
fpest <- coef_lst[[j]]$Pest[1]

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



fp_df$logit1min <- coef_lst[[j]]$intercept[1] + coef_lst[[j]]$slope[1] * log10(fp_df$minPD + 1)

fp_df$pred1min <- exp(fp_df$logit1min)/(1 + exp(fp_df$logit1min))

# if predicted value was too big, division by Inf has given NA - convert to 1s
fp_df$pred1min <- ifelse(is.na(fp_df$pred1min), 1, fp_df$pred1min)


fp_df$logit2min <- coef_lst[[j]]$intercept[2] + coef_lst[[j]]$slope[2] * log10(fp_df$minPD + 1)

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


# line types
mod2tips <- unique(coef_2a$FocalHost)
mod2notinc <- setdiff(phy$tip.label, mod2tips)


edgelines <- rep(1, Nedge(phy))

for(i in seq_along(mod2notinc)){
edgelines[which.edge(phy, mod2notinc[i])] <- 2
}



# plot ------------------------------------------------------------------------

#depth <- max(node.depth.edgelength(phy))

par(mfrow=c(1,2))
plot(phy,
     show.tip.label = FALSE,
#     x.lim = depth-(depth/9),
#     label.offset = 5,
     no.margin = TRUE,
     edge.color = edgecols1)
plot(phy,
     direction = "leftwards",
     no.margin = TRUE,
     show.tip.label = FALSE,
     edge.color = edgecols2,
     edge.lty = edgelines)




