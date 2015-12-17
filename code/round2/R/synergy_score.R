#!/usr/bin/env Rscript


## libraries
library(kernlab)


# Arguments
args=commandArgs(trailingOnly = TRUE)
#kernel_train_path=args[1]
#kernel_test_path=args[2]

#kernel_train=read.table(kernel_train_path,sep="\t",header=F)
#kernel_test=read.table(kernel_test_path,sep="\t",header=F)

############################################################ debugging
train_score=read.table("../../../data/originals/train_score.csv",sep=",",header=F)
kernel_train=read.table("../../../data/round2/kernel_train_test/kernel_train_dot_product_drug_pathway_angular_similarity_mutations_driver.txt",sep="\t",header=F)
kernel_test=read.table("../../../data/round2/kernel_train_test/kernel_test_dot_product_drug_pathway_angular_similarity_mutations_driver.txt",sep="\t",header=F)
###########################################################

########################################################### parameter range

############################################################ Cross Validation for ksvm #############################
num_train <- nrow(kernel_train)
k <- 10
num_fold <- round(num_train/k)
random_set <- sample(1:num_train, num_train, replace = FALSE)
training_index <- 1
validation_index <- -1
temp_num <- 0
validation_error_ksvr <- array(0,dim=k)
for(i in 1:k){
  index <- array(training_index,dim=num_train)
  if(i == k){
    index[random_set[(temp_num+1):num_train]] <- validation_index
  }else{
    index[random_set[(temp_num+1):(temp_num+num_fold)]] <- validation_index
  }
  kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[which(index==training_index),which(index==training_index)]))
  kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[which(index==validation_index),which(index==training_index)]))
  label_train <- as.matrix(train_score[which(index==training_index),1])
  label_validation <- as.matrix(train_score[which(index==validation_index),1])
  reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=1,type="nu-svr",cross=0,C = 50,nu = 1)
  sv_index <- reg_ksvm@alphaindex
  score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
  validation_error_ksvr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
  temp_num <- temp_num + num_fold
}
err <- sum(validation_error_ksvr)/k
print(err)

########################################################### cross validation for klr ##################################


validation_error_krr <- array(0,dim=k)
temp_num <- 0
for(i in 1:k){
  index <- array(training_index,dim=num_train)
  if(i == k){
    index[random_set[(temp_num+1):num_train]] <- validation_index
  }else{
    index[random_set[(temp_num+1):(temp_num+num_fold)]] <- validation_index
  }
  
  kernel_train_new <- as.matrix(kernel_train[which(index==training_index),which(index==training_index)])
  kernel_validation <- as.matrix(kernel_train[which(index==validation_index),which(index==training_index)])
  label_train <- as.matrix(train_score[which(index==training_index),1])
  label_validation <- as.matrix(train_score[which(index==validation_index),1])
  lambda <- 1000
  score_validation <- kernel_validation %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train 
  validation_error_krr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
  temp_num <- temp_num + num_fold
}
err <- sum(validation_error_krr)/k
print(err)




