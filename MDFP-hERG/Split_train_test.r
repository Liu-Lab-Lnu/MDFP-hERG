IC50 <- read.csv("IC50-8960-gai-2.csv", stringsAsFactors  = FALSE, header = TRUE)
IC50_data1 <- IC50[1:203,]
IC50_data2 <- IC50[204:nrow(IC50),]
dim(IC50)
dim(IC50_data1)
dim(IC50_data2)
nsample <- nrow(IC50_data2)
train_index <- sample(1:nsample,nsample*0.8)
IC50_data2_train <- IC50_data2[train_index,]
IC50_data2_test <- IC50_data2[-train_index,]
IC50_data2_train$Training.Test <- "Training"
IC50_data2_test$Training.Test <- "Test"

IC50_all <- rbind(IC50_data1, IC50_data2_train, IC50_data2_test)

dim(IC50_all)
summary(IC50)
summary(IC50_all)
table(IC50_all$Training.Test)

write.csv(IC50_all, "IC50_all.csv", row.names=FALSE)

