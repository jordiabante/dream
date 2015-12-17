#!/usr/bin/env Rscript


## libraries
library(kernlab)


# Arguments
args=commandArgs(trailingOnly = TRUE)
#kernel_train_path=args[1]
#kernel_test_path=args[2]
#kernel_test_path=args[3]

#kernel_train=read.table(kernel_train_path,sep="\t",header=F)
#kernel_test=read.table(kernel_test_path,sep="\t",header=F)

############################################################ debugging
train_score=read.table("../../../data/originals/train_score.csv",sep=",",header=F)
kernel_train=read.table("../../../data/round2/kernel_train_test/kernel_train_dot_product_drug_pathway_angular_similarity_mutations_driver.txt",sep="\t",header=F)
kernel_test=read.table("../../../data/round2/kernel_train_test/kernel_test_dot_product_drug_pathway_angular_similarity_mutations_driver.txt",sep="\t",header=F)
###########################################################
suffix <- strsplit(kernel_test_path,"./../../data/round2/kernel_train_test/kernel_test",)[[1]][2]
num_train <- nrow(kernel_train)
random_set <- sample(1:num_train, num_train, replace = FALSE)
num_validation <- 500
training_index <- 1
validation_index <- -1
temp_num <- 0
index <- array(training_index,dim=num_train)
index[random_set[1:num_validation]] <- validation_index
kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[which(index==training_index),which(index==training_index)]))
kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[which(index==validation_index),which(index==training_index)]))
label_train <- as.matrix(train_score[which(index==training_index),1])
label_validation <- as.matrix(train_score[which(index==validation_index),1])
reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=0.01,type="nu-svr",cross=0,C = 100,nu = 0.6)
sv_index <- reg_ksvm@alphaindex
score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
validation_error_ksvr <-sum((score_validation-label_validation)^2)/length(label_validation)
out1 <- cbind(paste("ksvr_",suffix,sep=""),t(reg_ksvm@fitted),t(score_validation))

########################################################### cross validation for klr ##################################

lambda <- 0.1
score_validation <- kernel_validation %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train
score_validation_train <- kernel_train_new %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train
validation_error_krr <-sum((score_validation-label_validation)^2)/length(label_validation)
out2 <- cbind(paste("ksvr_",suffix,sep=""),t(score_validation_train),t(score_validation))




