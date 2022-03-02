library(tidyverse)
library(caret)
source("scatterplot.r")

IC50 <- read.csv("IC50_all.csv", stringsAsFactors  = FALSE, header = TRUE)
pIC50 <- IC50


MDFP <- read.csv("MDFP_4_7836.csv", stringsAsFactors  = FALSE, header = TRUE)
MDFP3D11 <- read.csv("MDFP3D11-2-7837.csv", stringsAsFactors  = FALSE, header = TRUE)
BaselineFP <- read.csv("BaselineFP.csv", stringsAsFactors  = FALSE, header = TRUE)
FP3DFP <- read.csv("3D11-7837.csv", stringsAsFactors  = FALSE, header = TRUE)
PropertyFP <- read.csv("PropertyFP.csv", stringsAsFactors  = FALSE, header = TRUE)
ECFP4 <- read.csv("ECFP4.csv", stringsAsFactors  = FALSE, header = TRUE)

MDFP <- left_join(pIC50, MDFP, by = "Name")
MDFP3D11 <- left_join(pIC50, MDFP3D11, by = "Name")
BaselineFP <- left_join(pIC50, BaselineFP, by = "Name")
FP3DFP <- left_join(pIC50, FP3DFP, by = "Name")
PropertyFP <- left_join(pIC50, PropertyFP, by = "Name")
ECFP4 <- left_join(pIC50, ECFP4, by = "Name")

ALLFP <- left_join(MDFP, select(MDFP3D11, -pIC50, -Training.Test), by = "Name")
ALLFP <- left_join(ALLFP, select(BaselineFP, -pIC50, -Training.Test), by = "Name")
ALLFP <- left_join(ALLFP, select(FP3DFP, -pIC50, -Training.Test), by = "Name")
ALLFP <- left_join(ALLFP, select(PropertyFP, -pIC50, -Training.Test), by = "Name")
ALLFP <- left_join(ALLFP, select(ECFP4, -pIC50, -Training.Test), by = "Name")

MDFP_train <- filter(ALLFP, Training.Test == "Training") %>% select(-Name, -Training.Test)
MDFP_test <- filter(ALLFP, Training.Test == "Test") %>% select(-Name, -Training.Test)

#Choose features seleted by RFE
MDFP_rfe_predictors <- read.csv(file = "MDFP_rfe_predictors.csv", row.names = 1)
MDFP3D11_rfe_predictors <- read.csv(file = "MDFP3D11_rfe_predictors.csv", row.names = 1)
BaselineFP_rfe_predictors <- read.csv(file = "BaselineFP_rfe_predictors.csv", row.names = 1)
FP3DFP11_rfe_predictors <- read.csv(file = "3DFP11_rfe_predictors.csv", row.names = 1)
PropertyFP_rfe_predictors <- read.csv(file = "PropertyFP_rfe_predictors.csv", row.names = 1)
ECFP4_rfe_predictors <- read.csv(file = "ECFP4_rfe_predictors.csv", row.names = 1)

rfe_selected_predictors <- c(as.character(MDFP_rfe_predictors$x), as.character(MDFP3D11_rfe_predictors$x), as.character(BaselineFP_rfe_predictors$x), as.character(FP3DFP11_rfe_predictors$x), as.character(PropertyFP_rfe_predictors$x), as.character(ECFP4_rfe_predictors$x))
MDFP_train <- select(MDFP_train, pIC50, all_of(rfe_selected_predictors))


#Calculate the correlation of the selected features and pIC50
cor_pearson <- cor(MDFP_train, MDFP_train$pIC50, method = "pearson")
cor_kendall <- cor(MDFP_train,  MDFP_train$pIC50, method = "kendall")
cor_spearman <- cor(MDFP_train,  MDFP_train$pIC50, method = "spearman")

colnames(cor_pearson)[1] <- "pIC50"
colnames(cor_kendall)[1] <- "pIC50"
colnames(cor_spearman)[1] <- "pIC50"

write.csv(cor_pearson, file="cor_pearson.csv")
write.csv(cor_kendall, file="cor_kendall.csv")
write.csv(cor_spearman, file="cor_spearman.csv")

