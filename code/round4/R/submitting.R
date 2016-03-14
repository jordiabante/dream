#!/usr/bin/env Rscript


## libraries
library(kernlab)
library(quadprog)
# Arguments
#args=commandArgs(trailingOnly = TRUE)
############################################################# function
our_PCR <- function(A,y){
  PC <- princomp(A)
  K <- choosing_cut_off(A)
  beta <- solve(t(PC$scores[,1:K])%*% PC$scores[,1:K]+0.001 * diag(K))%*% t(PC$scores[,1:K]) %*% y
  alpha <- t(beta) %*% t(PC$loadings[,1:K])
  return(alpha)
}

choosing_cut_off <- function(A){
  num_train <- nrow(A)
  random_set <- sample(1:num_train, num_train, replace = FALSE)
  k <- 5
  num_fold <- round(num_train/k)
  training_index <- 1
  validation_index <- -1
  temp_num <- 0
  error <- matrix(0,nrow = 1,ncol = ncol(A))
  for(j in 1:k){
    index <- array(training_index,dim=num_train)
    if(j == k){
      index[random_set[(temp_num+1):num_train]] <- validation_index
    }else{
      index[random_set[(temp_num+1):(temp_num+num_fold)]] <- validation_index
    }
    train_new_set <- which(index==training_index)
    validation_set <- which(index==validation_index)
    train_score_new <- as.matrix(train_score[train_new_set,1])
    A_new <- A[train_new_set,]
    PC_new <- princomp(A_new)
    for(K in 1:ncol(A_new)){
      beta <- solve(t(PC_new$scores[,1:K])%*% PC_new$scores[,1:K]+0.001 * diag(K))%*% t(PC_new$scores[,1:K]) %*% train_score_new
      f_new <- t(beta) %*% t(PC_new$scores[,1:K])
      error[1,K] <- error[1,K]+ sum((t(f_new) - train_score_new)^2)
    }
    temp_num <- temp_num + num_fold
  }
  index <- which.min(error)
  return(index)
}


############################################################ debugging
parameter <- read.table("../../data/round4/param/param_global_product_non_weighted.txt",header = F,sep = "\t")
ch1_train <- read.csv("../../data/originals/ch1_train_combination_and_monoTherapy.csv", stringsAsFactors=F)
ch1_train <- ch1_train[which(ch1_train[,13]==1),]
ch1_test <- read.csv("../../data/originals/ch1_test_monoTherapy.csv", stringsAsFactors=F)
train_score <- as.data.frame(ch1_train[,12])
train_mean <- mean(train_score[,1])
train_sd <- sd(train_score[,1])
train_score <- as.matrix((train_score - train_mean)/(train_sd))
################## partintioning
random_set <- sample(1:nrow(train_score), nrow(train_score), replace = FALSE)
training_percent <- 0.7
num_train <- round(training_percent*nrow(train_score))
train_index <- random_set[1:num_train]
num_validation <- nrow(train_score)-num_train
validation_index <- random_set[(num_train+1):nrow(train_score)]
label_train <- train_score[train_index,1]
label_validation <- train_score[validation_index,1]
###################################

num_weak_learner <- nrow(parameter)
num_test <- nrow(ch1_test)
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,cbind(paste("ksvr_",parameter[i,1],sep=""),paste("krr_",parameter[i,1],sep="")))
}

weak_leaner_matrix_train <- matrix(0,nrow = num_train,ncol = 2*num_weak_learner)
weak_leaner_matrix_validation <- matrix(0,nrow = num_validation,ncol = 2*num_weak_learner)
weak_leaner_matrix_test <- matrix(0,nrow = num_test,ncol = 2*num_weak_learner)
colnames(weak_leaner_matrix_train) <- name_col
colnames(weak_leaner_matrix_test) <- name_col
colnames(weak_leaner_matrix_validation) <- name_col
for(i in 1:nrow(parameter)){
  print(i)
  kernel_train_original=read.table(paste("../../data/round4/kernel_train_test_product/kernel_train_",parameter[i,1],sep=""),sep="\t",header=F) 
  kernel_test=read.table(paste("../../data/round4/kernel_train_test_product/kernel_final_test_",parameter[i,1],sep=""),sep="\t",header=F)
  kernel_train <- as.kernelMatrix(as.matrix(kernel_train_original[train_index,train_index]))
  kernel_validation <- as.kernelMatrix(as.matrix(kernel_train_original[validation_index,train_index]))
  kernel_test <- as.kernelMatrix(as.matrix(kernel_test[,train_index]))
  reg_ksvm <- ksvm(kernel_train,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
  sv_index <- reg_ksvm@alphaindex
  score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
  score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
  score_test[which(is.na(score_test)==TRUE)] <- 0.25
  score_validation[which(is.na(score_validation)==TRUE)] <- 0.25
  weak_leaner_matrix_train[,2*i-1] <- reg_ksvm@fitted
  weak_leaner_matrix_test[,2*i-1] <- score_test
  weak_leaner_matrix_validation[,2*i-1] <- score_validation
  
  lambda <- parameter[i,5]
  ########
  #s <- vector(length = num_train)
  #index_less <- which(label_train < -2.25)
  #index_high <- which(label_train > 3.75)
  #index_other <-  which((label_train <= 3.75) & (label_train >= -2.25) )
  #s[index_less] <- 1
  #s[index_high] <- 14
#  s[index_other] <- ceiling(2*(label_train[index_other]+2.25))+1
  #my_hist <- hist(label_train, breaks = seq(from = -7.25,to = 8.25,by = 0.5))
  #count <- vector(length = 14)
  #count[1] <- sum(my_hist$counts[1:10])
  #count[14] <- sum(my_hist$counts[23:length(my_hist$counts)])
  #count[2:13] <- my_hist$counts[11:22]
  #d <- 1/count[s[1:length(s)]]
  #sqrt_D <- sqrt(length(s)) * diag(sqrt(d))/sqrt(sum(d))
  sqrt_D <- diag(length(label_train))
  ##########
  score_train_ridge <- kernel_train %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train)) + sqrt_D %*% kernel_train %*% sqrt_D) %*% sqrt_D %*% label_train
  score_test_ridge <- kernel_test %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train)) + sqrt_D %*% kernel_train %*% sqrt_D) %*% sqrt_D %*% label_train
  score_validation_ridge <- kernel_validation %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train)) + sqrt_D %*% kernel_train %*% sqrt_D) %*% sqrt_D %*% label_train
  score_test_ridge[which(is.na(score_test_ridge)==TRUE)] <- 0.25
  score_validation_ridge[which(is.na(score_validation_ridge)==TRUE)] <- 0.25
  weak_leaner_matrix_train[,2*i] <- score_train_ridge
  weak_leaner_matrix_test[,2*i] <- score_test_ridge
  weak_leaner_matrix_validation[,2*i] <- score_validation_ridge
}
################################################### QP
index_test <- unique(which(weak_leaner_matrix_test[,]>10,arr.ind = T)[,2])
index_train <- unique(which(weak_leaner_matrix_train[,]>10,arr.ind = T)[,2])
index_validation <- unique(which(weak_leaner_matrix_validation[,]>10,arr.ind = T)[,2])
bad_index <- unique(c(index_test,index_train,index_validation))

ensemble_train <- weak_leaner_matrix_validation[,-bad_index]
ensemble_label <- label_validation


R<-matrix(0,ncol(ensemble_train),ncol(ensemble_train))
for(l1 in 1:ncol(ensemble_train)){
  for(l2 in 1:ncol(ensemble_train)){
    R[l1,l2] <- sum((ensemble_train[,l1]-ensemble_label)%*%(ensemble_train[,l2]-ensemble_label))
  }
}
dvec <- matrix(0,1,ncol(ensemble_train))
I_mat       <- matrix(0,ncol(ensemble_train),ncol(ensemble_train))
diag(I_mat) <- 1
temp <- matrix(1,1,ncol(ensemble_train))
Amat <- t(rbind(temp[1,],I_mat))
bvec <- c(1,dvec)
Qp <- solve.QP(R,dvec,Amat,bvec,meq=1)

Qp_score_train <- t(Qp$solution %*% t(weak_leaner_matrix_train[,-bad_index]))
Qp_score_validation <- t(Qp$solution %*% t(weak_leaner_matrix_validation[,-bad_index]))
Qp_score_test <- t(Qp$solution %*% t(weak_leaner_matrix_test[,-bad_index]))
Qp_score_test <- Qp_score_test*train_sd+train_mean
result_name <- c("CELL_LINE","COMBINATION_ID","PREDICTION")

ch1_A <- cbind(ch1_test[,c(1,14)],Qp_score_test)

colnames(ch1_A) <- result_name
write.table(ch1_A,"../../data/final_submitting/ch1_A.csv",sep = ",",quote = F,row.names = F, col.names = T)
write.table(ch1_A,"../../data/final_submitting/ch1_A.txt",sep = "\t",quote = F,row.names = F, col.names = T)




mutation_index <- grep(pattern = "mutation",x = parameter[,1])
cnv_index <- grep(pattern = "cnv",x = parameter[,1])
index <- c(mutation_index,cnv_index)
part_2_index <- c()
for(i in index){
  part_2_index <- c(part_2_index,2*i-1,2*i)
}
weak_leaner_matrix_test <- weak_leaner_matrix_test[,part_2_index]
weak_leaner_matrix_train <- weak_leaner_matrix_train[,part_2_index]
weak_leaner_matrix_validation <- weak_leaner_matrix_validation[,part_2_index]

index_test <- unique(which(weak_leaner_matrix_test[,]>10,arr.ind = T)[,2])
index_train <- unique(which(weak_leaner_matrix_train[,]>10,arr.ind = T)[,2])
index_validation <- unique(which(weak_leaner_matrix_validation[,]>10,arr.ind = T)[,2])
bad_index <- unique(c(index_test,index_train,index_validation))

ensemble_train <- weak_leaner_matrix_validation[,-bad_index]
ensemble_label <- label_validation

R<-matrix(0,ncol(ensemble_train),ncol(ensemble_train))
for(l1 in 1:ncol(ensemble_train)){
  for(l2 in 1:ncol(ensemble_train)){
    R[l1,l2] <- sum((ensemble_train[,l1]-ensemble_label)%*%(ensemble_train[,l2]-ensemble_label))
  }
}
dvec <- matrix(0,1,ncol(ensemble_train))
I_mat       <- matrix(0,ncol(ensemble_train),ncol(ensemble_train))
diag(I_mat) <- 1
temp <- matrix(1,1,ncol(ensemble_train))
Amat <- t(rbind(temp[1,],I_mat))
bvec <- c(1,dvec)
Qp <- solve.QP(R,dvec,Amat,bvec,meq=1)

Qp_score_train <- t(Qp$solution %*% t(weak_leaner_matrix_train[,-bad_index]))
Qp_score_validation <- t(Qp$solution %*% t(weak_leaner_matrix_validation[,-bad_index]))
Qp_score_test <- t(Qp$solution %*% t(weak_leaner_matrix_test[,-bad_index]))
Qp_score_test <- Qp_score_test*train_sd+train_mean
result_name <- c("CELL_LINE","COMBINATION_ID","PREDICTION")

ch1_B <- cbind(ch1_test[,c(1,14)],Qp_score_test)

colnames(ch1_B) <- result_name
write.table(ch1_B,"../../data/final_submitting/ch1_B.csv",sep = ",",quote = F,row.names = F, col.names = T)
write.table(ch1_B,"../../data/final_submitting/ch1_B.txt",sep = "\t",quote = F,row.names = F, col.names = T)





