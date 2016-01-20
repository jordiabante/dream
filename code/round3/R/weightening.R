#library(neuralnet)
library(glm2)

num_validation <- 500
predictor_matrix <- read.table("../../../data/round2/predictor_matrix.txt",sep="\t",header = T)
response_predictor <- read.table("../../../data/round2/train_score_shuffle.txt",sep="\t",header = T)
predictor_matrix <- predictor_matrix[,1:100]
train_matrix <- predictor_matrix[(num_validation+1):nrow(predictor_matrix),]
test_matrix <- predictor_matrix[1:num_validation,]
train_label <- response_predictor[(num_validation+1):nrow(predictor_matrix),1]
test_label <- response_predictor[1:num_validation,1]
train <- cbind(train_matrix,train_label)
n <- colnames(train)[1:ncol(predictor_matrix)]
f <- as.formula(paste("train_label ~", paste(n[!n %in% "train_label"], collapse = " + ")))
########################

#nn <- neuralnet(f,data=train,hidden=c(3,2),linear.output=T)
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
glm_model$coefficients=glm_model$coefficients[!is.na(glm_model$coefficients)]
glm_predicted <- out$fitted.values
train_predicted_glm <- predict(glm_model,train)
train_error <-sum(abs(train_predicted_glm-train_label))/nrow(train)
test_predicted_glm <- predict(glm_model,test_matrix)
test_error <-sum(abs(test_predicted_glm-test_label))/nrow(test_matrix)
