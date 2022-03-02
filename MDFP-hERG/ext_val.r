load("../ALLFP-ML.rData")

library(tidyverse)
library(caret)
source("../scatterplot.r")

#Read data
IC50 <- read.csv("IC50-938.csv", stringsAsFactors  = FALSE, header = TRUE)
pIC50 <- IC50


MDFP <- read.csv("MDFP-1-938.csv", stringsAsFactors  = FALSE, header = TRUE)
MDFP3D11 <- read.csv("MDFP-2-938.csv", stringsAsFactors  = FALSE, header = TRUE)
BaselineFP <- read.csv("Baseline2D-938.csv", stringsAsFactors  = FALSE, header = TRUE)
FP3DFP <- read.csv("3DFP-938.csv", stringsAsFactors  = FALSE, header = TRUE)
PropertyFP <- read.csv("PropertyFP-938.csv", stringsAsFactors  = FALSE, header = TRUE)
ECFP4 <- read.csv("ECFP4-938.csv", stringsAsFactors  = FALSE, header = TRUE)

MDFP <- left_join(pIC50, MDFP, by = "Name")
MDFP3D11 <- left_join(pIC50, MDFP3D11, by = "Name")
BaselineFP <- left_join(pIC50, BaselineFP, by = "Name")
FP3DFP <- left_join(pIC50, FP3DFP, by = "Name")
PropertyFP <- left_join(pIC50, PropertyFP, by = "Name")
ECFP4 <- left_join(pIC50, ECFP4, by = "Name")

ALLFP <- left_join(MDFP, select(MDFP3D11, -pIC50), by = "Name")
ALLFP <- left_join(ALLFP, select(BaselineFP, -pIC50), by = "Name")
ALLFP <- left_join(ALLFP, select(FP3DFP, -pIC50), by = "Name")
ALLFP <- left_join(ALLFP, select(PropertyFP, -pIC50), by = "Name")
ALLFP <- left_join(ALLFP, select(ECFP4, -pIC50), by = "Name")


#Choose features seleted by RFE
MDFP_rfe_predictors <- read.csv(file = "../MDFP_rfe_predictors.csv", row.names = 1)
MDFP3D11_rfe_predictors <- read.csv(file = "../MDFP3D11_rfe_predictors.csv", row.names = 1)
BaselineFP_rfe_predictors <- read.csv(file = "../BaselineFP_rfe_predictors.csv", row.names = 1)
FP3DFP11_rfe_predictors <- read.csv(file = "../3DFP11_rfe_predictors.csv", row.names = 1)
PropertyFP_rfe_predictors <- read.csv(file = "../PropertyFP_rfe_predictors.csv", row.names = 1)
ECFP4_rfe_predictors <- read.csv(file = "../ECFP4_rfe_predictors.csv", row.names = 1)

rfe_selected_predictors <- c(as.character(MDFP_rfe_predictors$x), as.character(MDFP3D11_rfe_predictors$x), as.character(BaselineFP_rfe_predictors$x), as.character(FP3DFP11_rfe_predictors$x), as.character(PropertyFP_rfe_predictors$x), as.character(ECFP4_rfe_predictors$x))
ext_val_ALLFP <- select(ALLFP, pIC50, all_of(rfe_selected_predictors))

ext_val_pIC50 <- ext_val_ALLFP$pIC50


	#External Validation
	svmFit_pred <- predict(svmFit, newdata=ext_val_ALLFP)
	metric_test_svm <- postResample(pred = svmFit_pred, obs = ext_val_pIC50)
	plot_scatter(pred=svmFit_pred, obs=ext_val_pIC50, FP="ALLFP", ML=paste0("SVM", "_extALL"))


	gbmFit_pred <- predict(gbmFit, newdata=ext_val_ALLFP)
	metric_test_gbm <- postResample(pred = gbmFit_pred, obs = ext_val_pIC50)
	plot_scatter(pred=gbmFit_pred, obs=ext_val_pIC50, FP="ALLFP", ML=paste0("GBM", "_extALL"))


	rfFit_pred <- predict(rfFit, newdata=ext_val_ALLFP)
	metric_test_rf <- postResample(pred = rfFit_pred, obs = ext_val_pIC50)
	plot_scatter(pred=rfFit_pred, obs=ext_val_pIC50, FP="ALLFP", paste0("RF", "_extALL"))


	kknnFit_pred <- predict(kknnFit, newdata=ext_val_ALLFP)
	metric_test_kknn <- postResample(pred = kknnFit_pred, obs = ext_val_pIC50)
	plot_scatter(pred=kknnFit_pred, obs=ext_val_pIC50, FP="ALLFP", ML=paste0("kknn", "_extALL"))


	
	consensus_pred <- (svmFit_pred + gbmFit_pred + rfFit_pred + kknnFit_pred) /4
	metric_test_consensus <- postResample(pred = consensus_pred, obs = ext_val_pIC50)
	plot_scatter(pred=consensus_pred, obs=ext_val_pIC50, FP="ALLFP", ML=paste0("consensus", "_extALL"))

	consensus_pred_2 <- (svmFit_pred + gbmFit_pred + rfFit_pred) /3
	metric_test_consensus_2 <- postResample(pred = consensus_pred_2, obs = ext_val_pIC50)
	plot_scatter(pred=consensus_pred_2, obs=ext_val_pIC50, FP="ALLFP", ML=paste0("consensus_2", "_extALL"))

	#Output results
	metric_test <- bind_rows(svm=metric_test_svm, gbm=metric_test_gbm, rf=metric_test_rf, kknn=metric_test_kknn, consensus=metric_test_consensus, consensus_2=metric_test_consensus_2, .id = "method") %>% 
					select(method, RMSE, Rsquared, MAE)
					


write.csv(metric_test, file="ALLFP_metric_extALL.csv")

