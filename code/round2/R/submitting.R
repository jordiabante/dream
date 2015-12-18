#!/usr/bin/env Rscript


## libraries
library(kernlab)
library(glm2)

# Arguments
#args=commandArgs(trailingOnly = TRUE)


############################################################ debugging
label_train=read.table("../../../data/originals/train_score.csv",sep=",",header=F)
target_predictors <- read.table("../../../data/round2/target_predictors.txt",header=F)
pathway_predictors <- read.table("../../../data/round2/pathway_predictors.txt",header=F)


num_data  <- 2386
num_train <- 1795
parameter <-  target_predictors
synergy_score_matrix_target <- matrix(0,nrow = num_data,ncol = nrow(parameter))
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,paste("ksvr_",parameter[i,1],sep=""))
}
colnames(synergy_score_matrix_target) <- name_col
for(i in 1:nrow(parameter)){
  
  kernel_train=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_train_",parameter[i,1],".txt",sep=""),sep="\t",header=F))) 
  kernel_test=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_test_",parameter[i,1],".txt",sep=""),sep="\t",header=F)))
  reg_ksvm <- ksvm(kernel_train,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
  synergy_score_matrix_target[,i] <- rbind(score_test,reg_ksvm@fitted)
}


parameter <-  pathway_predictors
synergy_score_matrix_pathway <- matrix(0,nrow = num_data,ncol = nrow(parameter))
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,paste("ksvr_",parameter[i,1],sep=""))
}
colnames(synergy_score_matrix_pathway) <- name_col
for(i in 1:nrow(parameter)){
  kernel_train=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_train_",parameter[i,1],".txt",sep=""),sep="\t",header=F))) 
  kernel_test=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_test_",parameter[i,1],".txt",sep=""),sep="\t",header=F)))
  reg_ksvm <- ksvm(kernel_train,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
  synergy_score_matrix_pathway[,i] <- rbind(score_test,reg_ksvm@fitted)
}
################################################################### first submission: target
train_matrix <- synergy_score_matrix_target[1:num_train,]
test_matrix <- synergy_score_matrix_target[(num_train+1):num_data,]
train <- cbind(train_matrix,label_train)
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_target_a <- predict(glm_model,as.data.frame(test_matrix))
##########################################################################
train_matrix <- synergy_score_matrix_pathway[1:num_train,]
test_matrix <- synergy_score_matrix_pathway[(num_train+1):num_data,]
train <- cbind(train_matrix,label_train)
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_pathway_a <- predict(glm_model,as.data.frame(test_matrix))
##########################################################################
train_matrix <- synergy_score_matrix_pathway[1:num_train,]
test_matrix <- synergy_score_matrix_pathway[(num_train+1):num_data,]
train <- cbind(synergy_score_matrix_pathway[1:num_train,],synergy_score_matrix_target[1:num_train,],label_train)
test <- cbind(synergy_score_matrix_pathway[(num_train+1):num_data,],synergy_score_matrix_target[(num_train+1):num_data,])
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_all_a <- predict(glm_model,as.data.frame(test))





target_predictors <- read.table("../../../data/round2/target_predictors_b.txt",header=F)
pathway_predictors <- read.table("../../../data/round2/pathway_predictors_b.txt",header=F)


parameter <-  target_predictors
synergy_score_matrix_target <- matrix(0,nrow = num_data,ncol = nrow(parameter))
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,paste("ksvr_",parameter[i,1],sep=""))
}
colnames(synergy_score_matrix_target) <- name_col
for(i in 1:nrow(parameter)){
  
  kernel_train=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_train_",parameter[i,1],".txt",sep=""),sep="\t",header=F))) 
  kernel_test=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_test_",parameter[i,1],".txt",sep=""),sep="\t",header=F)))
  reg_ksvm <- ksvm(kernel_train,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
  synergy_score_matrix_target[,i] <- rbind(score_test,reg_ksvm@fitted)
}


parameter <-  pathway_predictors
synergy_score_matrix_pathway <- matrix(0,nrow = num_data,ncol = nrow(parameter))
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,paste("ksvr_",parameter[i,1],sep=""))
}
colnames(synergy_score_matrix_pathway) <- name_col
for(i in 1:nrow(parameter)){
  kernel_train=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_train_",parameter[i,1],".txt",sep=""),sep="\t",header=F))) 
  kernel_test=as.kernelMatrix(as.matrix(read.table(paste("../../../data/round2/kernel_train_test/kernel_test_",parameter[i,1],".txt",sep=""),sep="\t",header=F)))
  reg_ksvm <- ksvm(kernel_train,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
  synergy_score_matrix_pathway[,i] <- rbind(score_test,reg_ksvm@fitted)
}
################################################################### first submission: target
train_matrix <- synergy_score_matrix_target[1:num_train,]
test_matrix <- synergy_score_matrix_target[(num_train+1):num_data,]
train <- cbind(train_matrix,label_train)
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_target_b <- predict(glm_model,as.data.frame(test_matrix))
##########################################################################
train_matrix <- synergy_score_matrix_pathway[1:num_train,]
test_matrix <- synergy_score_matrix_pathway[(num_train+1):num_data,]
train <- cbind(train_matrix,label_train)
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_pathway_b <- predict(glm_model,as.data.frame(test_matrix))
##########################################################################
train_matrix <- synergy_score_matrix_pathway[1:num_train,]
test_matrix <- synergy_score_matrix_pathway[(num_train+1):num_data,]
train <- cbind(synergy_score_matrix_pathway[1:num_train,],synergy_score_matrix_target[1:num_train,],label_train)
test <- cbind(synergy_score_matrix_pathway[(num_train+1):num_data,],synergy_score_matrix_target[(num_train+1):num_data,])
n <- colnames(train)[1:ncol(synergy_score_matrix_target)]
f <- as.formula(paste("V1 ~", paste(n[!n %in% "V1"], collapse = " + ")))
glm_model <- glm(f, data = train,family = gaussian(link = "identity"))
predicted_glm_all_b <- predict(glm_model,as.data.frame(test))















