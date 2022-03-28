# LJ 2022-03-09 plot model predictions presentably on phylogeny



# setup -----------------------------------------------------------------------

library(ape)
library(viridis)


# generate tip-level predictions for plotting ---------------------------------

for(j in 1:length(coef_lst)){


# choose a pest
#j=1

fpest <- coef_lst[[j]]$Pest[1]

# extract focal pest column from big db
fp_df <- data.frame(t(db[fpest,]))

colnames(fp_df) <- "compatible"

# subset to only potential hosts in phy
#fp_df <- fp_df[which(rownames(fp_df)%in%phy$tip.label),]

# id known hosts
hosts <- rownames(fp_df)[fp_df$compatible==1]

if(length(hosts)>1){

# get PD all tips->known hosts
for(i in seq_along(hosts)){

tmp <- cophenetic(phy)[,hosts[i]]

fp_df[,ncol(fp_df)+1] <- tmp[match(rownames(fp_df), names(tmp))]

colnames(fp_df)[ncol(fp_df)] <- paste0("PDfh", i)
}

# average PD
fp_df$avPD <- rowMeans(fp_df[,-which(colnames(fp_df)=="compatible")])

# shortest PD to any known host
fp_df$minPD <- apply(fp_df[,-which(colnames(fp_df)%in%c("compatible","avPD"))],
                     1, min)


# use first set of coefficients to get predicted values for all tips based on
# min PD->known host
fp_df$logit1min <- coef_lst[[j]]$intercept[1] + coef_lst[[j]]$slope[1] * log10(fp_df$minPD + 1)

# covert to p(compatible)
fp_df$pred1min <- exp(fp_df$logit1min)/(1 + exp(fp_df$logit1min))

# if predicted value was too big, division by Inf has given NA - convert to 1s
fp_df$pred1min <- ifelse(is.na(fp_df$pred1min), 1, fp_df$pred1min)


# repeat for 2nd set of coefficients
fp_df$logit2min <- coef_lst[[j]]$intercept[2] + coef_lst[[j]]$slope[2] * log10(fp_df$minPD + 1)

fp_df$pred2min <- exp(fp_df$logit2min)/(1 + exp(fp_df$logit2min))

fp_df$pred2min <- ifelse(is.na(fp_df$pred2min), 1, fp_df$pred2min)



# plotting palette
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
edgelines <- rep(1, Nedge(phy))

misTips <- colnames(db_2)[which(is.na(db_2[fpest,]))]

for(i in seq_along(misTips)){
edgelines[which.edge(phy, misTips[i])] <- 2
}



# plot ------------------------------------------------------------------------

phydir <- paste0("output/phyplots/", runID, "_vs", subID)

if(!dir.exists(phydir)){
dir.create(phydir, recursive = TRUE)
}

png(paste0(phydir, "/pest", j, ".png"),
           width = 400, height = 600, unit = "px")


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

dev.off()
}else{print(paste0("only one host for pest ", fpest, " (", j, ")"))}

}

