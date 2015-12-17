#!/usr/bin/env Rscript


## libraries
library(kernlab)


# Arguments
#args=commandArgs(trailingOnly = TRUE)


############################################################ debugging
train_score=read.table("../../../data/originals/train_score.csv",sep=",",header=F)
parameter <- read.table("../../../data/round2/cross_validation_optimal.txt",header=F)
num_train <- 1795
random_set <- sample(1:num_train, num_train, replace = FALSE)
num_validation <- 500
training_index <- 1
validation_index <- -1
index <- array(training_index,dim=num_train)
index[random_set[1:num_validation]] <- validation_index
label_train <- as.matrix(train_score[which(index==training_index),1])
label_validation <- as.matrix(train_score[which(index==validation_index),1])
synergy_score_matrix <- matrix(0,nrow = num_train,ncol = 2*nrow(parameter))
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,cbind(paste("ksvr_",parameter[i,1],sep=""),paste("krr_",parameter[i,1],sep="")))
}
colnames(synergy_score_matrix) <- name_col
for(i in 1:nrow(parameter)){
  print(i)
  kernel_train=read.table(paste("../../../data/round2/kernel_train_test/kernel_train_",parameter[i,1],".txt",sep=""),sep="\t",header=F) 
  kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[which(index==training_index),which(index==training_index)]))
  kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[which(index==validation_index),which(index==training_index)]))
  reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
  validation_error_ksvr <-sum((score_validation-label_validation)^2)/length(label_validation)
  synergy_score_matrix[,2*i-1] <- rbind(score_validation,reg_ksvm@fitted)
  
  
  lambda <- parameter[i,6]
  score_validation <- kernel_validation %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train
  score_validation_train <- kernel_train_new %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train
  validation_error_krr <-sum((score_validation-label_validation)^2)/length(label_validation)
  synergy_score_matrix[,2*i] <- rbind(score_validation,score_validation_train)
  
}
write.table(synergy_score_matrix,"../../../data/round2/predictor_matrix.txt",sep="\t",col.names = T,row.names = F,quote = F)



