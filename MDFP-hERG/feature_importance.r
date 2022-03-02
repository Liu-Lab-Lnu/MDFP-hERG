library(tidyverse)
library(randomForest)

IC50 <- read.csv("IC50_all.csv", stringsAsFactors  = FALSE, header = TRUE)
pIC50 <- IC50

MDFP <- read.csv("MDFP_4_7836.csv", stringsAsFactors  = FALSE, header = TRUE)
PropertyFP <- read.csv("PropertyFP.csv", stringsAsFactors  = FALSE, header = TRUE)
BaselineFP <- read.csv("BaselineFP.csv", stringsAsFactors  = FALSE, header = TRUE)
ECFP4 <- read.csv("ECFP4.csv", stringsAsFactors  = FALSE, header = TRUE)
MDFP3D11 <- read.csv("MDFP3D11-2-7837.csv", stringsAsFactors  = FALSE, header = TRUE)
FP3DFP <- read.csv("3D11-7837.csv", stringsAsFactors  = FALSE, header = TRUE)


MDFP <- left_join(pIC50, MDFP, by = "Name")
PropertyFP <- left_join(pIC50, PropertyFP, by = "Name")
BaselineFP <- left_join(pIC50, BaselineFP, by = "Name")
ECFP4 <- left_join(pIC50, ECFP4, by = "Name")
MDFP3D11 <- left_join(pIC50, MDFP3D11, by = "Name")
FP3DFP <- left_join(pIC50, FP3DFP, by = "Name")


MDFP_train <- filter(MDFP, Training.Test == "Training") %>% select(-Name, -Training.Test)
PropertyFP_train <- filter(PropertyFP, Training.Test == "Training") %>% select(-Name, -Training.Test)
BaselineFP_train <- filter(BaselineFP, Training.Test == "Training") %>% select(-Name, -Training.Test)
ECFP4_train <- filter(ECFP4, Training.Test == "Training") %>% select(-Name, -Training.Test)
MDFP3D11_train <- filter(MDFP3D11, Training.Test == "Training") %>% select(-Name, -Training.Test)
FP3DFP_train <- filter(FP3DFP, Training.Test == "Training") %>% select(-Name, -Training.Test)

#Choose features seleted by RFE
MDFP_rfe_predictors <- read.csv(file = "MDFP_rfe_predictors.csv", row.names = 1)
PropertyFP_rfe_predictors <- read.csv(file = "PropertyFP_rfe_predictors.csv", row.names = 1)
BaselineFP_rfe_predictors <- read.csv(file = "BaselineFP_rfe_predictors.csv", row.names = 1)
ECFP4_rfe_predictors <- read.csv(file = "ECFP4_rfe_predictors.csv", row.names = 1)
MDFP3D11_rfe_predictors <- read.csv(file = "MDFP3D11_rfe_predictors.csv", row.names = 1)
FP3DFP_rfe_predictors <- read.csv(file = "3DFP11_rfe_predictors.csv", row.names = 1)

MDFP_train <- select(MDFP_train, pIC50, all_of(as.character(MDFP_rfe_predictors$x)))
PropertyFP_train <- select(PropertyFP_train, pIC50, all_of(as.character(PropertyFP_rfe_predictors$x)))
BaselineFP_train <- select(BaselineFP_train, pIC50, all_of(as.character(BaselineFP_rfe_predictors$x)))
ECFP4_train <- select(ECFP4_train, pIC50, all_of(as.character(ECFP4_rfe_predictors$x)))
MDFP3D11_train <- select(MDFP3D11_train, pIC50, all_of(as.character(MDFP3D11_rfe_predictors$x)))
FP3DFP_train <- select(FP3DFP_train, pIC50, all_of(as.character(FP3DFP_rfe_predictors$x)))

#Use Random Forest to calculate the importance of the selected features
MDFP_rf <- randomForest(pIC50 ~ ., data = MDFP_train, importance = TRUE, ntree=500)
PropertyFP_rf <- randomForest(pIC50 ~ ., data = PropertyFP_train, importance = TRUE)
BaselineFP_rf <- randomForest(pIC50 ~ ., data = BaselineFP_train, importance = TRUE)
ECFP4_rf <- randomForest(pIC50 ~ ., data = ECFP4_train, importance = TRUE)
MDFP3D11_rf <- randomForest(pIC50 ~ ., data = MDFP3D11_train, importance = TRUE)
FP3DFP_rf <- randomForest(pIC50 ~ ., data = FP3DFP_train, importance = TRUE)

#Ptint the results of importance 
write.csv(importance(MDFP_rf, scale=TRUE) , file = "MDFP_feature_importance.csv")
write.csv(importance(PropertyFP_rf, scale=TRUE) , file = "PropertyFP_feature_importance.csv")
write.csv(importance(BaselineFP_rf, scale=TRUE), file = "BaselineFP_feature_importance.csv")
write.csv(importance(ECFP4_rf, scale=TRUE) , file = "ECFP4_feature_importance.csv")
write.csv(importance(MDFP3D11_rf, scale=TRUE) , file = "MDFP3D11_feature_importance.csv")
write.csv(importance(FP3DFP_rf, scale=TRUE) , file = "FP3DFP_feature_importance.csv")


#Output significance plot results
pdf(file = "MDFP_feature_importance.pdf")
varImpPlot(MDFP_rf)
dev.off()

pdf(file = "PropertyFP_feature_importance.pdf")
varImpPlot(PropertyFP_rf)
dev.off()


pdf(file = "BaselineFP_feature_importance.pdf")
varImpPlot(BaselineFP_rf)
dev.off()


pdf(file = "ECFP4_feature_importance.pdf")
varImpPlot(ECFP4_rf)
dev.off()

pdf(file = "MDFP3D11_feature_importance.pdf")
varImpPlot(MDFP3D11_rf)
dev.off()

pdf(file = "FP3DFP_feature_importance.pdf")
varImpPlot(FP3DFP_rf)
dev.off()






