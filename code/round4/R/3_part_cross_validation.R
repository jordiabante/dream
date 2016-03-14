#!/usr/bin/env Rscript


## libraries
library(kernlab)
library(quadprog)
library(Matrix)
library(rdetools)
library(plyr)
# Arguments
#args=commandArgs(trailingOnly = TRUE)
############################################################# function
inner_angular_product <- function(feature){
  my_kernel_inner <- diag(rep(1,nrow(feature)))
  my_kernel_angular <- diag(rep(1,nrow(feature)))
  for (i in 1:(nrow(feature)-1)){
    for (j in (i+1):nrow(feature)){
      my_kernel_inner[i,j] <- (feature[i,] %*% feature[j,])/(norm(feature[i,],"2")*norm(feature[j,],"2"))
      my_kernel_angular[i,j] <- 1-acos(pmin(pmax(my_kernel_inner[i,j],-1.0),1.0))/pi
      #my_kernel_angular[i,j] <- 1-2*acos(pmin(pmax(my_kernel_inner[i,j],-1.0),1.0))/pi
    }
  }
  my_kernel_inner <- as.matrix(forceSymmetric(as.matrix(my_kernel_inner)))
  my_kernel_angular <- as.matrix(forceSymmetric(as.matrix(my_kernel_angular)))
  return(rbind(my_kernel_inner,my_kernel_angular))
}

global_error <- function(x,y,obs){
  agg <- aggregate(SYNERGY_SCORE ~ CELL_LINE, obs, median)
  z0 <- agg$SYNERGY_SCORE[match(obs$CELL_LINE, agg$CELL_LINE)]
  
  agg <- aggregate(SYNERGY_SCORE ~ COMBINATION_ID, obs, median)
  z1 <- agg$SYNERGY_SCORE[match(obs$COMBINATION_ID, agg$COMBINATION_ID)]  
  numerator=parCor(x,y,z1) - parCor(x,z0,z1) * parCor(z0,y,z1)
  denumerator= sqrt(1-parCor(x,z0,z1)^2) * sqrt(1-parCor(z0,y,z1)^2)
  return(c(numerator/denumerator))
}

final_error <- function(x,y,obs){
  agg <- aggregate(SYNERGY_SCORE ~ CELL_LINE, obs, median)
  z0 <- agg$SYNERGY_SCORE[match(obs$CELL_LINE, agg$CELL_LINE)]
  
  agg <- aggregate(SYNERGY_SCORE ~ COMBINATION_ID, obs, median)
  z1 <- agg$SYNERGY_SCORE[match(obs$COMBINATION_ID, agg$COMBINATION_ID)]  
  numerator=parCor(x,y,z1) - parCor(x,z0,z1) * parCor(z0,y,z1)
  denumerator= sqrt(1-parCor(x,z0,z1)^2) * sqrt(1-parCor(z0,y,z1)^2)
  ##############
  temp   <- data.frame(OBSERVATION = x, PREDICTION = y, COMBINATION_ID = obs$COMBINATION_ID)
  R      <- ddply(temp, "COMBINATION_ID", function(x){if(length(x$OBSERVATION) > 1){return(cor(x$OBSERVATION,x$PREDICTION))}else{return(0)}})
  R[is.na(R[,2]),2] <- 0;
  R$N    <- table(obs$COMBINATION_ID)
  obsMax <- ddply(obs, "COMBINATION_ID", function(x){max(x$SYNERGY_SCORE)})
  
  primaryScoreCH1     <- sum(R$V1*sqrt(R$N-1))/sum(sqrt(R$N-1))                       # primary metric
  gt20Inds             <- obsMax[,2] >= 20
  tieBreakingScoreCH1  <- sum((R$V1*sqrt(R$N-1))[gt20Inds])/sum(sqrt(R$N-1)[gt20Inds]) # tie-breaking metric
  # ----------------------------------------
  
  # partial out the mean of synergy across cell lines and combinationations
  return(c(score=numerator/denumerator,
           final= primaryScoreCH1,
           tiebreak= tieBreakingScoreCH1))
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

our_PCR <- function(A,y){
  PC <- princomp(A)
  K <- choosing_cut_off(A)
  beta <- solve(t(PC$scores[,1:K])%*% PC$scores[,1:K])%*% t(PC$scores[,1:K]) %*% y
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
    label_train <- as.matrix(train_score[train_new_set,1])
    label_validation <- as.matrix(train_score[validation_set,1])
    A_new <- A[train_new_set,]
    PC_new <- princomp(A_new)
    for(K in 1:ncol(A_new)){
      beta <- solve(t(PC_new$scores[,1:K])%*% PC_new$scores[,1:K]+0.001 * diag(K))%*% t(PC_new$scores[,1:K]) %*% label_train
      #f_new <- t(beta) %*% t(PC_new$scores[,1:K])
      #error[1,K] <- error[1,K]+ sum((t(f_new) - label_train)^2)
      alpha <- t(beta) %*% t(PC_new$loadings[,1:K])
      f_test <- t(alpha %*% t(A[validation_set,]))
      error[1,K] <- error[1,K]+ sum((f_test - label_validation)^2)
    }
    temp_num <- temp_num + num_fold
  }
  index <- which.min(error)
  return(index)
}

############################################################ should be changed
parameter <- read.table("../../data/round4/param/param_global_product_non_weighted.txt",header=F,sep = "\t")
ch1_train <- read.csv("../../data/originals/ch1_train_combination_and_monoTherapy.csv", stringsAsFactors=F)
##################################
#parameter <- parameter[1:18,]
ch1_train <- ch1_train[which(ch1_train[,13]==1),]
train_score <- as.data.frame(ch1_train[,12])
our_mean <-  mean(train_score[,1])
our_sd <- sd(train_score[,1])
train_score <- (train_score - our_mean)/(our_sd)
num_train <- nrow(train_score)
random_set <- sample(1:num_train, num_train, replace = FALSE)
k <- 4
num_fold <- round(num_train/k)
num_weak_learner <- 2*nrow(parameter)
training_index <- 1
validation_index <- -1
test_index <- 2
num_method <- 3
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,cbind(paste("ksvr_",parameter[i,1],sep=""),paste("krr_",parameter[i,1],sep="")))
}
cross_validation_error <- array(0,dim = c(12,k*(k-1),num_weak_learner+num_method))
temp_num <- 0
num_cv <- 0
for(j in 1:(k-1)){
  for(l in 1:(2*(k-j))){
    index <- array(training_index,dim=num_train)
    if(l %% 2 == 0){
      index[random_set[(temp_num+1):(temp_num+num_fold)]] <- test_index
      if( l == 2*(k-j)){
        index[random_set[(temp_num+ceiling(l/2)*num_fold+1):(num_train)]] <- validation_index
      }else{
        index[random_set[(temp_num+ceiling(l/2)*num_fold+1):(temp_num+(ceiling(l/2)+1)*num_fold)]] <- validation_index
      }
    }else{
      index[random_set[(temp_num+1):(temp_num+num_fold)]] <- validation_index
      if( l == 2*(k-j)-1){
        index[random_set[(temp_num+ceiling(l/2)*num_fold+1):(num_train)]] <- test_index
      }else{
        index[random_set[(temp_num+ceiling(l/2)*num_fold+1):(temp_num+(ceiling(l/2)+1)*num_fold)]] <- test_index
      }
    }
    num_cv <- num_cv +1
    train_new_set <- which(index==training_index)
    validation_set <- which(index==validation_index)
    test_set <- which(index==test_index)
    label_train <- as.matrix(train_score[train_new_set,1])
    label_validation <- as.matrix(train_score[validation_set,1])
    label_test <- as.matrix(train_score[test_set,1])
    weak_leaner_matrix_train <- matrix(0,nrow = length(train_new_set),ncol = num_weak_learner)
    weak_leaner_matrix_validation <- matrix(0,nrow = length(validation_set),ncol = num_weak_learner)
    weak_leaner_matrix_test <- matrix(0,nrow = length(test_set),ncol = num_weak_learner)
    colnames(weak_leaner_matrix_train) <- name_col
    colnames(weak_leaner_matrix_validation) <- name_col
    colnames(weak_leaner_matrix_test) <- name_col
    for(i in 1:nrow(parameter)){
      print(i) 
      kernel_train=read.table(paste("../../data/round4/kernel_train_test_product/kernel_train_",parameter[i,1],sep=""),sep="\t",header=F)       
      #############################
      kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[train_new_set,train_new_set]))
      kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[validation_set,train_new_set]))
      kernel_test <- as.kernelMatrix(as.matrix(kernel_train[test_set,train_new_set]))
      reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
      sv_index <- reg_ksvm@alphaindex
      score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
      score_test<- predict(reg_ksvm,as.kernelMatrix(kernel_test[,sv_index]))
      temp_train <- final_error(label_train[,1],as.numeric(reg_ksvm@fitted),ch1_train[train_new_set,c(1,14,12)])
      temp_validation <-  final_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
      temp_test <-  final_error(label_test[,1],as.numeric(score_test),ch1_train[test_set,c(1,14,12)])
      cross_validation_error[1,num_cv,2*i-1] <- reg_ksvm@error
      cross_validation_error[2,num_cv,2*i-1] <-sum((score_validation-label_validation)^2)/length(label_validation)
      cross_validation_error[3,num_cv,2*i-1] <- sum((score_test-label_test)^2)/length(label_test)
      cross_validation_error[4,num_cv,2*i-1] <- temp_train[1]
      cross_validation_error[5,num_cv,2*i-1] <- temp_validation[1]
      cross_validation_error[6,num_cv,2*i-1] <- temp_test[1]
      cross_validation_error[7,num_cv,2*i-1] <- temp_train[2]
      cross_validation_error[8,num_cv,2*i-1] <- temp_validation[2]
      cross_validation_error[9,num_cv,2*i-1] <- temp_test[2]
      cross_validation_error[10,num_cv,2*i-1] <- temp_train[3]
      cross_validation_error[11,num_cv,2*i-1] <- temp_validation[3]
      cross_validation_error[12,num_cv,2*i-1] <- temp_test[3]
      
      weak_leaner_matrix_train[,2*i-1] <- reg_ksvm@fitted
      weak_leaner_matrix_validation[,2*i-1] <- score_validation
      weak_leaner_matrix_test[,2*i-1] <- score_test
      
      lambda <- parameter[i,5]
      ####################################### first version
      #s <- vector(length = nrow(label_train))
      #index_less <- which(label_train < -2.25)
      #index_high <- which(label_train > 3.75)
      #index_other <-  which((label_train <= 3.75) & (label_train >= -2.25) )
      #s[index_less] <- 1
      #s[index_high] <- 14
      #s[index_other] <- ceiling(2*(label_train[index_other]+2.25))+1
      #my_hist <- hist(label_train, breaks = seq(from = -7.25,to = 8.25,by = 0.5))
      #count <- vector(length = 14)
      #count[1] <- sum(my_hist$counts[1:10])
      #count[14] <- sum(my_hist$counts[23:length(my_hist$counts)])
      #count[2:13] <- my_hist$counts[11:22]
      #d <- 1/count[s[1:length(s)]]
      #sqrt_D <- sqrt(length(s)) * diag(sqrt(d))/sqrt(sum(d))
      ###################################### second version
      #s <- vector(length = nrow(label_train))
      #index_less <- which(label_train < -0.25)
      #index_high <- which(label_train > 0.25)
      #index_other <-  which((label_train <= .25) & (label_train >= -0.25) )
      #s[index_less] <- 1
      #s[index_high] <- 2
      #s[index_other] <- 3
      #my_hist <- hist(label_train, breaks = c(-7.25,0.25,-0.25,8.25))
      #count <- my_hist$counts
      #d <- 1/count[s[1:length(s)]]
      #sqrt_D <- sqrt(length(s)) * diag(sqrt(d))/sqrt(sum(d))
      sqrt_D <- diag(nrow(label_train))
      
      score_validation <- kernel_validation %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
      score_validation_train <- kernel_train_new %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
      score_test <- kernel_test %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
      
      temp_train <- final_error(label_train[,1],as.numeric(score_validation_train),ch1_train[train_new_set,c(1,14,12)])
      temp_validation <-  final_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
      temp_test <-  final_error(label_test[,1],as.numeric(score_test),ch1_train[test_set,c(1,14,12)])
      
      cross_validation_error[1,num_cv,2*i] <-sum((score_validation_train-label_train)^2)/length(label_train)
      cross_validation_error[2,num_cv,2*i] <-sum((score_validation-label_validation)^2)/length(label_validation)
      cross_validation_error[3,num_cv,2*i] <-sum((score_test-label_test)^2)/length(label_test)
      cross_validation_error[4,num_cv,2*i] <- temp_train[1]
      cross_validation_error[5,num_cv,2*i] <- temp_validation[1]
      cross_validation_error[6,num_cv,2*i] <- temp_test[1]
      cross_validation_error[7,num_cv,2*i] <- temp_train[2]
      cross_validation_error[8,num_cv,2*i] <- temp_validation[2]
      cross_validation_error[9,num_cv,2*i] <- temp_test[2]
      cross_validation_error[10,num_cv,2*i] <- temp_train[3]
      cross_validation_error[11,num_cv,2*i] <- temp_validation[3]
      cross_validation_error[12,num_cv,2*i] <- temp_test[3]      
      
      weak_leaner_matrix_train[,2*i] <- score_validation_train
      weak_leaner_matrix_validation[,2*i] <- score_validation
      weak_leaner_matrix_test[,2*i] <- score_test
    }
    ##################################################################################### not seen data
    ############################### glm
    method_counter <- 0
    ensemble_train <- weak_leaner_matrix_validation
    ensemble_label <- label_validation
    #ensemble_train_df <- as.data.frame(cbind(ensemble_train,ensemble_label))
    #f <- as.formula(paste("ensemble_label ~", paste(name_col[!name_col %in% "ensemble_label"], collapse = " + ")))
    #glm_model <- glm(f, data = ensemble_train_df,family = gaussian(link = "identity"))
    #train_predicted_glm <- predict(glm_model,as.data.frame(weak_leaner_matrix_train))
    #validation_predicted_glm <- predict(glm_model,as.data.frame(weak_leaner_matrix_validation))
    #test_predicted_glm <- predict(glm_model,as.data.frame(weak_leaner_matrix_test))
    
    #temp_train <- final_error(label_train[,1],as.numeric(train_predicted_glm),ch1_train[train_new_set,c(1,14,12)])
    #temp_validation <-  final_error(label_validation[,1],as.numeric(validation_predicted_glm),ch1_train[validation_set,c(1,14,12)])
    #temp_test <-  final_error(label_test[,1],as.numeric(test_predicted_glm),ch1_train[test_set,c(1,14,12)])
    #
    #cross_validation_error[1,num_cv,2*i+method_counter] <-sum(abs(train_predicted_glm-label_train)^2)/nrow(label_train)
    #cross_validation_error[2,num_cv,2*i+method_counter] <-sum(abs(validation_predicted_glm-label_validation)^2)/nrow(weak_leaner_matrix_validation)
    #cross_validation_error[3,num_cv,2*i+method_counter] <-sum(abs(test_predicted_glm-label_test)^2)/nrow(weak_leaner_matrix_test)
    #cross_validation_error[4,num_cv,2*i+method_counter] <- temp_train[1]
    #cross_validation_error[5,num_cv,2*i+method_counter] <- temp_validation[1]
    #cross_validation_error[6,num_cv,2*i+method_counter] <- temp_test[1]
    #cross_validation_error[7,num_cv,2*i+method_counter] <- temp_train[2]
    #cross_validation_error[8,num_cv,2*i+method_counter] <- temp_validation[2]
    #cross_validation_error[9,num_cv,2*i+method_counter] <- temp_test[2]
    #cross_validation_error[10,num_cv,2*i+method_counter] <- temp_train[3]
    #cross_validation_error[11,num_cv,2*i+method_counter] <- temp_validation[3]
    #cross_validation_error[12,num_cv,2*i+method_counter] <- temp_test[3]    
    
    
    ############################### PCR
    method_counter <- method_counter+1
    PCR_weight <- our_PCR(as.matrix(ensemble_train),ensemble_label)
    PCR_score_train <- t(PCR_weight %*% t(weak_leaner_matrix_train))
    PCR_score_validation <- t(PCR_weight %*% t(weak_leaner_matrix_validation))
    PCR_score_test <- t(PCR_weight %*% t(weak_leaner_matrix_test))
    
    temp_train <- final_error(label_train[,1],as.numeric(PCR_score_train),ch1_train[train_new_set,c(1,14,12)])
    temp_validation <-  final_error(label_validation[,1],as.numeric(PCR_score_validation),ch1_train[validation_set,c(1,14,12)])
    temp_test <-  final_error(label_test[,1],as.numeric(PCR_score_test),ch1_train[test_set,c(1,14,12)])
    
    
    cross_validation_error[1,num_cv,2*i+method_counter]<-sum(abs(as.matrix(PCR_score_train)-label_train)^2)/nrow(weak_leaner_matrix_train)
    cross_validation_error[2,num_cv,2*i+method_counter] <-sum(abs(as.matrix(PCR_score_validation)-label_validation)^2)/nrow(weak_leaner_matrix_validation)
    cross_validation_error[3,num_cv,2*i+method_counter] <-sum(abs(as.matrix(PCR_score_test)-label_test)^2)/nrow(weak_leaner_matrix_test)
    cross_validation_error[4,num_cv,2*i+method_counter] <- temp_train[1]
    cross_validation_error[5,num_cv,2*i+method_counter] <- temp_validation[1]
    cross_validation_error[6,num_cv,2*i+method_counter] <- temp_test[1]
    cross_validation_error[7,num_cv,2*i+method_counter] <- temp_train[2]
    cross_validation_error[8,num_cv,2*i+method_counter] <- temp_validation[2]
    cross_validation_error[9,num_cv,2*i+method_counter] <- temp_test[2]
    cross_validation_error[10,num_cv,2*i+method_counter] <- temp_train[3]
    cross_validation_error[11,num_cv,2*i+method_counter] <- temp_validation[3]
    cross_validation_error[12,num_cv,2*i+method_counter] <- temp_test[3]    
    ############################### ridge
    method_counter <- method_counter+1
    ens_lamda <- 1
    ridge_score_validation <- weak_leaner_matrix_validation %*% solve(t(ensemble_train) %*% ensemble_train+ ens_lamda * diag(ncol(ensemble_train))) %*% t(ensemble_train) %*% ensemble_label
    ridge_score_train <- weak_leaner_matrix_train %*% solve(t(ensemble_train) %*% ensemble_train+ ens_lamda * diag(ncol(ensemble_train))) %*% t(ensemble_train) %*% ensemble_label
    ridge_score_test <- weak_leaner_matrix_test %*% solve(t(ensemble_train) %*% ensemble_train+ ens_lamda * diag(ncol(ensemble_train))) %*% t(ensemble_train) %*% ensemble_label
    
    temp_train <- final_error(label_train[,1],as.numeric(ridge_score_train),ch1_train[train_new_set,c(1,14,12)])
    temp_validation <-  final_error(label_validation[,1],as.numeric(ridge_score_validation),ch1_train[validation_set,c(1,14,12)])
    temp_test <-  final_error(label_test[,1],as.numeric(ridge_score_test),ch1_train[test_set,c(1,14,12)])
    
    cross_validation_error[1,num_cv,2*i+method_counter]<-sum(abs(ridge_score_train-label_train)^2)/nrow(weak_leaner_matrix_train)
    cross_validation_error[2,num_cv,2*i+method_counter] <-sum(abs(ridge_score_validation-label_validation)^2)/nrow(weak_leaner_matrix_validation)
    cross_validation_error[3,num_cv,2*i+method_counter] <-sum(abs(ridge_score_test-label_test)^2)/nrow(weak_leaner_matrix_test)
    cross_validation_error[4,num_cv,2*i+method_counter] <- temp_train[1]
    cross_validation_error[5,num_cv,2*i+method_counter] <- temp_validation[1]
    cross_validation_error[6,num_cv,2*i+method_counter] <- temp_test[1]
    cross_validation_error[7,num_cv,2*i+method_counter] <- temp_train[2]
    cross_validation_error[8,num_cv,2*i+method_counter] <- temp_validation[2]
    cross_validation_error[9,num_cv,2*i+method_counter] <- temp_test[2]
    cross_validation_error[10,num_cv,2*i+method_counter] <- temp_train[3]
    cross_validation_error[11,num_cv,2*i+method_counter] <- temp_validation[3]
    cross_validation_error[12,num_cv,2*i+method_counter] <- temp_test[3]
    ################################## Quadratic programming
    method_counter <- method_counter+1
    R<-matrix(0,2*nrow(parameter),2*nrow(parameter))
    for(l1 in 1:(2*nrow(parameter))){
      for(l2 in 1:(2*nrow(parameter))){
        R[l1,l2] <- sum((ensemble_train[,l1]-ensemble_label[,1])%*%(ensemble_train[,l2]-ensemble_label[,1]))
      }
    }
    dvec <- matrix(0,1,2*nrow(parameter))
    I_mat       <- matrix(0,2*nrow(parameter),2*nrow(parameter))
    diag(I_mat) <- 1
    temp <- matrix(1,1,2*nrow(parameter))
    Amat <- t(rbind(temp[1,],I_mat))
    bvec <- c(1,dvec)
    Qp <- solve.QP(R,dvec,Amat,bvec,meq=0)
    
    Qp_score_train <- t(Qp$solution %*% t(weak_leaner_matrix_train))
    Qp_score_validation <- t(Qp$solution %*% t(weak_leaner_matrix_validation))
    Qp_score_test <- t(Qp$solution %*% t(weak_leaner_matrix_test))
    
    temp_train <- final_error(label_train[,1],as.numeric(Qp_score_train),ch1_train[train_new_set,c(1,14,12)])
    temp_validation <-  final_error(label_validation[,1],as.numeric(Qp_score_validation),ch1_train[validation_set,c(1,14,12)])
    temp_test <-  final_error(label_test[,1],as.numeric(Qp_score_test),ch1_train[test_set,c(1,14,12)])
    
    cross_validation_error[1,num_cv,2*i+method_counter]<-sum(abs(Qp_score_train-label_train)^2)/nrow(weak_leaner_matrix_train)
    cross_validation_error[2,num_cv,2*i+method_counter] <-sum(abs(Qp_score_validation-label_validation)^2)/nrow(weak_leaner_matrix_validation)
    cross_validation_error[3,num_cv,2*i+method_counter] <-sum(abs(Qp_score_test-label_test)^2)/nrow(weak_leaner_matrix_test)
    cross_validation_error[4,num_cv,2*i+method_counter] <- temp_train[1]
    cross_validation_error[5,num_cv,2*i+method_counter] <- temp_validation[1]
    cross_validation_error[6,num_cv,2*i+method_counter] <- temp_test[1]
    cross_validation_error[7,num_cv,2*i+method_counter] <- temp_train[2]
    cross_validation_error[8,num_cv,2*i+method_counter] <- temp_validation[2]
    cross_validation_error[9,num_cv,2*i+method_counter] <- temp_test[2]
    cross_validation_error[10,num_cv,2*i+method_counter] <- temp_train[3]
    cross_validation_error[11,num_cv,2*i+method_counter] <- temp_validation[3]
    cross_validation_error[12,num_cv,2*i+method_counter] <- temp_test[3]
  }
  temp_num <- temp_num + num_fold
}
x <- cross_validation_error[9,,]
x <- apply(x, MARGIN=c(2), sum)/num_cv
plot(x)
#write.table(cross_validation_error,"cross_validation_error_high_dim.txt",quote = F,row.names = F,col.names = F,sep = "\t")
