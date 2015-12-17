#!/usr/bin/env Rscript


## libraries
library(kernlab)


# Arguments
args=commandArgs(trailingOnly = TRUE)
kernel_train_path=args[1]
kernel_test_path=args[2]

kernel_train=read.table(kernel_train_path,sep="\t",header=F)
kernel_test=read.table(kernel_test_path,sep="\t",header=F)

############################################################ debugging
train_score=read.table("../../data/originals/train_score.csv",sep=",",header=F)
#kernel_train=read.table("kernel_train_mutation.csv",sep=",",header=F)
#kernel_test=read.table("kernel_test_mutation.csv",sep=",",header=F)
###########################################################

########################################################### parameter range
eps_range <- c(0.01,0.1,0.2,0.5,0.7,1,2,3,5,10)
lambda_range <- c(0.005,0.05,0.1,0.5,1,10,30,50,70,100,500,1000) 
C_range <- c(0.005,0.05,0.1,0.5,1,10,30,50,70,100,500,1000)
nu_range <- c(0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
############################################################ Cross Validation for ksvm #############################
num_train <- nrow(kernel_train)
k <- 10
num_fold <- round(num_train/k)
random_set <- sample(1:num_train, num_train, replace = FALSE)
training_index <- 1
validation_index <- -1
best_max <- 10000000
best_index_eps <- 0
for(l in 1: length(eps_range)){
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
    reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=eps_range[l],type="nu-svr",cross=0,C = 50,nu = 1)
    sv_index <- reg_ksvm@alphaindex
    score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
    validation_error_ksvr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    temp_num <- temp_num + num_fold
  }
  err <- sum(validation_error_ksvr)/k
  if(err < best_max){
    bes_max <- err
    best_index_eps <- l
  }
}

best_max <- 10000000
best_index_C <- 0
for(l in 1: length(C_range)){
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
    reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=eps_range[best_index_eps],type="nu-svr",cross=0,C = C_range[l],nu = 1)
    sv_index <- reg_ksvm@alphaindex
    score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
    validation_error_ksvr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    temp_num <- temp_num + num_fold
  }
  err <- sum(validation_error_ksvr)/k
  if(err < best_max){
    bes_max <- err
    best_index_C <- l
  }
}

best_max <- 10000000
best_index_nu <- 0
for(l in 1: length(nu_range)){
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
    reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=eps_range[best_index_eps],type="nu-svr",cross=0,C = C_range[best_index_C],nu = nu_range[l])
    sv_index <- reg_ksvm@alphaindex
    score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
    validation_error_ksvr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    temp_num <- temp_num + num_fold
  }
  err <- sum(validation_error_ksvr)/k
  if(err < best_max){
    bes_max <- err
    best_index_nu <- l
  }
}

C_best <- C_range[best_index_C]
eps_best <- eps_range[best_index_eps]
nu_best <- nu_range[best_index_nu]
ksvr_err <- bes_max
########################################################### cross validation for klr ##################################

best_max <- 10000000
best_index_lamda <- 0
for(l in 1: length(lambda_range)){
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
    lambda <- lambda_range[l]
    score_validation <- kernel_validation %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train 
    validation_error_krr[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    temp_num <- temp_num + num_fold
  }
  err <- sum(validation_error_krr)/k
  if(err < best_max){
    bes_max <- err
    best_index_lambda <- l
  }
}
lamda_best <- lambda_range[best_index_lambda]
krrr_err <- bes_max
output <- cbind(kernel_train_path,eps_best,C_best,nu_best,round(ksvr_err),lamda_best,round(krrr_err))
write.table(output,stdout(),sep="\t",row.names=F,col.names=F,quote=F)


