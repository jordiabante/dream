#!/usr/bin/env Rscript

## libraries
library(kernlab)
library(quadprog)
library(Matrix)
library(rdetools)
library(plyr)

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
################################################################
ch1_train <- read.csv("../../data/originals/ch1_train_combination_and_monoTherapy.csv", stringsAsFactors=F,sep=",",header=T)
ch1_train <- ch1_train[which(ch1_train[,13]==1),]
train_score <- as.data.frame(ch1_train[,12])
train_score <- (train_score - mean(train_score[,1]))/(sd(train_score[,1]))
########################################################### parameter range
eps_range <- c(0.001,0.005,0.01,0.05,0.1,0.5,1,2,5,10,20)
C_range <- c(0.01,0.1,1,10,20,30,40,50,60,70,80,90,100,150,200,500,1000)
nu_range <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
#############################################################

num_train <- nrow(kernel_train)
k <- 5
num_fold <- round(num_train/k)
random_set <- sample(1:num_train, num_train, replace = FALSE)
training_index <- 1
validation_index <- -1
for(l1 in 1: length(eps_range)){
  for(l2 in 1: length(C_range)){
    for(l3 in 1: length(nu_range)){
      temp_num <- 0
      validation_error_ksvr_cor <- array(0,dim=k)
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
        reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=eps_range[l1],type="nu-svr",cross=0,C = C_range[l2],nu = nu_range[l3])
        sv_index <- reg_ksvm@alphaindex
        score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
        temp <-final_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
        validation_error_ksvr_cor[i] <- temp[2]
        temp_num <- temp_num + num_fold
      }
      err_cor <- sum(validation_error_ksvr_cor)/k
      output <- cbind(eps_range[l1],C_range[l2],nu_range[l3],err_cor)
      write.table(output,stdout(),sep="\t",row.names=F,col.names=F,quote=F)
    }
  }
}
