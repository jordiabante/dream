#!/usr/bin/env Rscript

library(kernlab)

# Arguments
args=commandArgs(trailingOnly = TRUE)
kernel_train_path=args[1]
kernel_test_path=args[2]

kernel_train=read.table(kernel_train_path,sep="\t",header=F)
kernel_test=read.table(kernel_test_path,sep="\t",header=F)
############################################################# function
global_error <- function(x,y,obs){
  agg <- aggregate(SYNERGY_SCORE ~ CELL_LINE, obs, median)
  z0 <- agg$SYNERGY_SCORE[match(obs$CELL_LINE, agg$CELL_LINE)]
  
  agg <- aggregate(SYNERGY_SCORE ~ COMBINATION_ID, obs, median)
  z1 <- agg$SYNERGY_SCORE[match(obs$COMBINATION_ID, agg$COMBINATION_ID)]  
  numerator=parCor(x,y,z1) - parCor(x,z0,z1) * parCor(z0,y,z1)
  denumerator= sqrt(1-parCor(x,z0,z1)^2) * sqrt(1-parCor(z0,y,z1)^2)
  return(c(numerator/denumerator))
}


parCor <- function(u,v,w) {
  numerator = cor(u,v) - cor(u,w) * cor(w,v)
  denumerator = sqrt(1-cor(u,w)^2) * sqrt(1-cor(w,v)^2)
  return(numerator/denumerator)
}

our_mean_error <- function(pre,obs){
  pred <- obs
  pred$SYNERGY_SCORE <- pre
  pred <- pred[match(paste(obs$CELL_LINE,obs$COMBINATION_ID),paste(pred$CELL_LINE,pred$COMBINATION_ID)),]
  R <- c()
  pred$COMBINATION_ID <- gsub(" ", "", pred$COMBINATION_ID)
  for (i in as.character(unique(obs$COMBINATION_ID))) {
    R <- c(R, cor(obs[obs$COMBINATION_ID == i, 'SYNERGY_SCORE'], 
                  pred[pred$COMBINATION_ID == i, 'SYNERGY_SCORE']))
  }
  #Make NA's in R = 0
  R[is.na(R)] = 0
  names(R) <- as.character(unique(obs$COMBINATION_ID))
  return(mean(R))
}
################################################################
ch1_train <- read.csv("../../data/originals/ch1_train_combination_and_monoTherapy.csv", stringsAsFactors=F,sep=",",header=T)
ch1_train <- ch1_train[which(ch1_train[,13]==1),]
train_score <- as.data.frame(ch1_train[,12])
train_score <- (train_score - mean(train_score[,1]))/(sd(train_score[,1]))
########################################################### parameter range
lambda_range <-  c(0.01,0.1,0.2,0.5,0.8,1,2,5,10,20,30,40,50,60,70,80,90,100,150,200,500,1000)
##############################################################
num_train <- nrow(kernel_train)
k <- 5
num_fold <- round(num_train/k)
random_set <- sample(1:num_train, num_train, replace = FALSE)
training_index <- 1
validation_index <- -1
# ksvr_err_mse <- best_max_mse
# ksvr_err_global <- best_max_global
# best_max_global <- 0
# best_max_mse <- 10000000
# index_lambda_global <- 0
# index_lambda_mse <- 0
# best_max_global_weight <- 0
# best_max_mse_weight <- 10000000
# index_lambda_global_weight <- 0
# index_lambda_mse_weight <- 0
for(l in 1: length(lambda_range)){
  validation_error_krr_mse <- array(0,dim=k)
  validation_error_krr_global <- array(0,dim=k)
  validation_error_krr_mean <- array(0,dim=k)
  validation_error_krr_mse_weight <- array(0,dim=k)
  validation_error_krr_global_weight <- array(0,dim=k)
  validation_error_krr_mean_weight <- array(0,dim=k)
  temp_num <- 0
  for(i in 1:k){
    index <- array(training_index,dim=num_train)
    if(i == k){
      index[random_set[(temp_num+1):num_train]] <- validation_index
    }else{
      index[random_set[(temp_num+1):(temp_num+num_fold)]] <- validation_index
    }
    train_new_set <- which(index==training_index)
    validation_set <- which(index==validation_index)
    label_train <- as.matrix(train_score[train_new_set,1])
    label_validation <- as.matrix(train_score[validation_set,1])
    kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[train_new_set,train_new_set]))
    kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[validation_set,train_new_set]))
    lambda <- lambda_range[l]
    ###############################################
    s <- vector(length = nrow(label_train))
    index_less <- which(label_train < -2.25)
    index_high <- which(label_train > 3.75)
    index_other <-  which((label_train <= 3.75) & (label_train >= -2.25) )
    s[index_less] <- 1
    s[index_high] <- 14
    s[index_other] <- ceiling(2*(label_train[index_other]+2.25))+1
    my_hist <- hist(label_train, breaks = seq(from = -7.25,to = 8.25,by = 0.5))
    count <- vector(length = 14)
    count[1] <- sum(my_hist$counts[1:10])
    count[14] <- sum(my_hist$counts[23:length(my_hist$counts)])
    count[2:13] <- my_hist$counts[11:22]
    d <- 1/count[s[1:length(s)]]
    sqrt_D <- sqrt(length(s)) * diag(sqrt(d))/sqrt(sum(d))
    score_validation_weight <- kernel_validation %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
    score_validation <- kernel_validation %*% solve(lambda * diag(nrow(kernel_train_new)) + kernel_train_new) %*% label_train 
    validation_error_krr_mse[i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    validation_error_krr_global[i] <-global_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
    validation_error_krr_mean[i] <-our_mean_error(score_validation,ch1_train[validation_set,c(1,14,12)])
    validation_error_krr_mse_weight[i] <-sum((score_validation_weight-label_validation)^2)/length(label_validation)
    validation_error_krr_global_weight[i] <-global_error(label_validation[,1],as.numeric(score_validation_weight),ch1_train[validation_set,c(1,14,12)])
    validation_error_krr_mean_weight[i] <-our_mean_error(score_validation_weight,ch1_train[validation_set,c(1,14,12)])
    temp_num <- temp_num + num_fold
  }
  err_mse <- sum(validation_error_krr_mse)/k
  #   if(err_mse < best_max_mse){
  #     print(err_mse)
  #     best_max_mse <- err_mse
  #     index_lambda_mse <- l
  #   }
  
  err_global <- sum(validation_error_krr_global)/k
  #   if(err_global > best_max_global){
  #     print(err_global)
  #     best_max_global <- err_global
  #     index_lambda_global <- l
  #   }
  
  err_mean <- sum(validation_error_krr_mean)/k
  
  err_mse_weight <- sum(validation_error_krr_mse_weight)/k
  #   if(err_mse_weight < best_max_mse_weight){
  #     print(err_mse_weight)
  #     best_max_mse_weight <- err_mse_weight
  #     index_lambda_mse_weight <- l
  #   }
  
  err_global_weight <- sum(validation_error_krr_global_weight)/k
  #   if(err_global_weight > best_max_global_weight){
  #     print(err_global_weight)
  #     best_max_global_weight <- err_global_weight
  #     index_lambda_global_weight <- l
  #   }
  err_mean_weight <- sum(validation_error_krr_mean_weight)/k
  output <- cbind(lambda_range[l],err_mse,err_global,err_mean,err_mse_weight,err_global_weight,err_mean_weight)
  write.table(output,stdout(),sep="\t",row.names=F,col.names=F,quote=F)
}
#}
#}
