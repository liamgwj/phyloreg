# LJ 2022-03-24 average model coefficients per pest group

coef <- read.csv("/home/liam/Documents/MSc/analysis/phyloreg/datasets/Parker/model_coefficients/cleaned_Parker_database_2022-02-28_p067_coef.csv")

coef$PestGroup <- gsub("[0-9]", "", coef$Pest)

coef_av <- aggregate(coef[,-which(colnames(coef)%in%c("Pest", "FocalHost",
                                                      "PestGroup"))],
                     list(Pest = coef$PestGroup),
                     mean)

write.csv(coef_av, "/home/liam/Documents/MSc/analysis/phyloreg/datasets/Parker/model_coefficients/cleaned_Parker_database_2022-02-28_p067_GroupCoefs.csv", row.names=FALSE)
