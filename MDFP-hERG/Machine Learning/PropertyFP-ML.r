library(tidyverse)
library(caret)
source("scatterplot.r")
#Read data
IC50 <- read.csv("IC50_all.csv", stringsAsFactors  = FALSE, header = TRUE)
pIC50 <- IC50


MDFP <- read.csv("PropertyFP.csv", stringsAsFactors  = FALSE, header = TRUE)

MDFP <- left_join(pIC50, MDFP, by = "Name")

MDFP_train <- filter(MDFP, Training.Test == "Training") %>% select(-Name, -Training.Test)
MDFP_test <- filter(MDFP, Training.Test == "Test") %>% select(-Name, -Training.Test)

# Remove Zero- and Near Zero-Variance Predictors
MDFP_train <- MDFP_train[,-(nearZeroVar(MDFP_train))]

#Multi-core parallel computing
library(doMC)
registerDoMC(cores = 20)
do_rfe_select <- FALSE
if (do_rfe_select){
	#RFE feature selection
	rfRFE <-  list(summary = defaultSummary,
				   fit = function(x, y, first, last, ...){
					 library(randomForest)
					 randomForest(x, y, importance = first, ...)
					 },
				   pred = function(object, x)  predict(object, x),
				   rank = function(object, x, y) {
					 vimp <- varImp(object)
					 vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
					 vimp$var <- rownames(vimp)                  
					 vimp
					 },
				   selectSize = pickSizeBest,
				   selectVar = pickVars)
	ctrl <- rfeControl(functions = rfRFE,
					   method = "repeatedcv",
					   number = 5,
					   repeats = 10,
					   verbose = TRUE)
	subsets <- seq(1, ncol(MDFP_train)-1)
	set.seed(2)
	rfeProfile <- rfe(pIC50 ~ ., data = MDFP_train,
					 sizes = subsets,
					 rfeControl = ctrl)

	write.csv(predictors(rfeProfile), file = "PropertyFP_rfe_predictors.csv")

	#Choose features seleted by RFE
	MDFP_train <- select(MDFP_train, pIC50, predictors(rfeProfile))
}else{
	PropertyFP_rfe_predictors <- read.csv(file = "PropertyFP_rfe_predictors.csv", row.names = 1)
	MDFP_train <- select(MDFP_train, pIC50, all_of(as.character(PropertyFP_rfe_predictors$x)))
}







#Use 10 different random seeds for training, 2^seq(1,10)
metric_cv_all <- data.frame()
metric_test_all <- data.frame()
for (seed in 2^seq(1,10)){
	print("seed:")
	print(seed)
	print("\n")


	#Five fold cross validation for 10 repetitions
	fitControl1 <- trainControl(method = "repeatedcv", 
							   number = 5, 
							   repeats = 1, 
							   search = 'random', 
							   verboseIter = TRUE)

	fitControl2 <- trainControl(method = "repeatedcv", 
							   number = 5, 
							   repeats = 1, 
							   search = 'grid', 
							   verboseIter = TRUE)

	svmGrid <- expand.grid(sigma = exp(1)^seq(-7,-3,1), 
						   C = exp(1)^seq(0,6,1))
	set.seed(seed)
	svmFit  <- train(pIC50 ~ ., data = MDFP_train, 
								method = "svmRadial", 
								preProc = c("center", "scale"), 
								tuneGrid = svmGrid,
								trControl = fitControl2,
								verbose = FALSE)
								
	gbmGrid <- expand.grid(n.trees = (1:100)*100, 
						   interaction.depth = 15, 
						   shrinkage = 0.02,
						   n.minobsinnode = 10)
	set.seed(seed)
	gbmFit  <- train(pIC50 ~ ., data = MDFP_train, 
								method = "gbm", 
								preProc = c("center", "scale"), 
								tuneGrid = gbmGrid,
								trControl = fitControl2,
								verbose = FALSE)

	rfGrid <- expand.grid(mtry = (max(floor((ncol(MDFP_train)-1)/3)-10,1)):min(floor((ncol(MDFP_train)-1)/3)+15, ncol(MDFP_train)-1))
	set.seed(seed)
	rfFit  <- train(pIC50 ~ ., data = MDFP_train, 
								method = "rf", 
								preProc = c("center", "scale"), 
								tuneGrid = rfGrid,
								ntree = 500,
								trControl = fitControl2,
								verbose = FALSE)

	kknnGrid <- expand.grid(kmax  = 30, distance=c(0.5, 1, 1.5, 2), kernel=c("rectangular", "triangular", "epanechnikov", "biweight", "triweight", "cos", "inv", "gaussian"))
	set.seed(seed)
	kknnFit  <- train(pIC50 ~ ., data = MDFP_train, 
								method = "kknn", 
								preProc = c("center", "scale"), 
								tuneGrid = kknnGrid,
								# tuneLength = 40,
								trControl = fitControl2)

	kknnGrid

	#Test set validation
	metric_cv_svm <- svmFit$results %>% top_n(1, desc(RMSE))
	svmFit_pred <- predict(svmFit, newdata=MDFP_test)
	metric_test_svm <- postResample(pred = svmFit_pred, obs = MDFP_test$pIC50)
	plot_scatter(pred=svmFit_pred, obs=MDFP_test$pIC50, FP="PropertyFP", ML=paste0("SVM", "_seed", seed))

	svm_grid_plot <- ggplot(svmFit$results, aes(log(sigma), log(C), fill=RMSE)) + geom_raster() + scale_fill_continuous(low="green", high="red")
	ggsave(svm_grid_plot, file = "tuneplots/svm_grid_plot_PropertyFP.pdf")


	metric_cv_gbm <- gbmFit$results %>% top_n(1, desc(RMSE))
	gbmFit_pred <- predict(gbmFit, newdata=MDFP_test)
	metric_test_gbm <- postResample(pred = gbmFit_pred, obs = MDFP_test$pIC50)
	plot_scatter(pred=gbmFit_pred, obs=MDFP_test$pIC50, FP="PropertyFP", ML=paste0("GBM", "_seed", seed))

	gbm_grid_plot <- ggplot(gbmFit$results, aes(n.trees, RMSE, color = as.factor(interaction.depth))) + geom_point() + geom_line()
	ggsave(gbm_grid_plot, file = "tuneplots/gbm_grid_plot_PropertyFP.pdf")



	metric_cv_rf <- rfFit$results %>% top_n(1, desc(RMSE))
	rfFit_pred <- predict(rfFit, newdata=MDFP_test)
	metric_test_rf <- postResample(pred = rfFit_pred, obs = MDFP_test$pIC50)
	plot_scatter(pred=rfFit_pred, obs=MDFP_test$pIC50, FP="PropertyFP", paste0("RF", "_seed", seed))

	rf_grid_plot <- ggplot(rfFit$results, aes(mtry, RMSE)) + geom_point() + geom_line()
	ggsave(rf_grid_plot, file = "tuneplots/rf_grid_plot_PropertyFP.pdf")


	metric_cv_kknn <- kknnFit$results %>% top_n(1, desc(RMSE))
	kknnFit_pred <- predict(kknnFit, newdata=MDFP_test)
	metric_test_kknn <- postResample(pred = kknnFit_pred, obs = MDFP_test$pIC50)
	plot_scatter(pred=kknnFit_pred, obs=MDFP_test$pIC50, FP="PropertyFP", ML=paste0("kknn", "_seed", seed))

	kknn_grid_plot <- ggplot(kknnFit$results, aes(distance, RMSE, color = kernel)) + geom_point() + geom_line() 
	ggsave(kknn_grid_plot, file = "tuneplots/kknn_grid_plot_PropertyFP.pdf")
	
	consensus_pred <- (svmFit_pred + gbmFit_pred + rfFit_pred + kknnFit_pred) /4
	metric_test_consensus <- postResample(pred = consensus_pred, obs = MDFP_test$pIC50)
	plot_scatter(pred=consensus_pred, obs=MDFP_test$pIC50, FP="PropertyFP", ML=paste0("consensus", "_seed", seed))

	consensus_pred_2 <- (svmFit_pred + gbmFit_pred + rfFit_pred) /3
	metric_test_consensus_2 <- postResample(pred = consensus_pred_2, obs = MDFP_test$pIC50)
	plot_scatter(pred=consensus_pred_2, obs=MDFP_test$pIC50, FP="PropertyFP", ML=paste0("consensus_2", "_seed", seed))

	#Combie results
	metric_cv <- bind_rows(svm=metric_cv_svm, gbm=metric_cv_gbm, rf=metric_cv_rf, kknn=metric_cv_kknn, .id = "method") %>%
					select(method, RMSE, Rsquared, MAE, RMSESD, RsquaredSD, MAESD) %>%
					mutate(seed=seed)
	metric_test <- bind_rows(svm=metric_test_svm, gbm=metric_test_gbm, rf=metric_test_rf, kknn=metric_test_kknn, consensus=metric_test_consensus, consensus_2=metric_test_consensus_2, .id = "method") %>% 
					select(method, RMSE, Rsquared, MAE) %>%
					mutate(seed=seed)
					
	metric_cv_all <- bind_rows(metric_cv_all, metric_cv)
	metric_test_all <- bind_rows(metric_test_all, metric_test)
}

#Summerized results
metric_cv_all_summary <- filter(metric_cv_all, RMSE < 2.0) %>%
							group_by(method) %>%
							summarise(across(RMSE:MAE, .fns=list(mean=mean,sd=sd)))

metric_test_all_summary <- filter(metric_test_all, RMSE < 2.0) %>%
							group_by(method) %>%
							summarise(across(RMSE:MAE, .fns=list(mean=mean,sd=sd)))


write.csv(metric_cv_all, file="PropertyFP_metric_cv_all.csv")
write.csv(metric_test_all, file="PropertyFP_metric_test_all.csv")

write.csv(metric_cv_all_summary, file="PropertyFP_metric_cv_all_summary.csv")
write.csv(metric_test_all_summary, file="PropertyFP_metric_test_all_summary.csv")

save.image("PropertyFP-ML.rData")

print("finished.")
