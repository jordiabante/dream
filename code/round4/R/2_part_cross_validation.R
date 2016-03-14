#!/usr/bin/env Rscript


## libraries
library(kernlab)
library(plyr)
#library(deepnet)
# Arguments
#args=commandArgs(trailingOnly = TRUE)
############################################################# function
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
parameter <- read.table("../../data/round4/param/param_global_product.txt",header=F,sep = "\t")

ch1_train <- read.csv("../../data/originals/ch1_train_combination_and_monoTherapy.csv", stringsAsFactors=F)
############################################################
ch1_train <- ch1_train[which(ch1_train[,13]==1),]
train_score <- as.data.frame(ch1_train[,12])
train_score <- (train_score - mean(train_score[,1]))/(sd(train_score[,1]))
num_train <- nrow(train_score)
random_set <- sample(1:num_train, num_train, replace = FALSE)
k <- 5
num_fold <- round(num_train/k)
num_weak_learner <- 2*nrow(parameter)
training_index <- 1
validation_index <- -1
num_method <- 3
temp_num <- 0
name_col <- c()
for(i in 1:nrow(parameter)){
  name_col <- cbind(name_col,cbind(paste("ksvr_",parameter[i,1],sep=""),paste("krr_",parameter[i,1],sep="")))
}
cross_validation_error <- array(0,dim = c(10,k,num_weak_learner+num_method))

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
  temp_num <- temp_num + num_fold
  weak_leaner_matrix_train <- matrix(0,nrow = length(train_new_set),ncol = num_weak_learner)
  weak_leaner_matrix_validation <- matrix(0,nrow = length(validation_set),ncol = num_weak_learner)
  colnames(weak_leaner_matrix_train) <- name_col
  colnames(weak_leaner_matrix_validation) <- name_col
  for(i in 1:nrow(parameter)){
    #if (i == 50){
    print(i)
    kernel_train=read.table(paste("../../dream/data/round3/kernel_train_test_product/kernel_train_",parameter[i,1],sep=""),sep="\t",header=F) 
    kernel_train_new <- as.kernelMatrix(as.matrix(kernel_train[train_new_set,train_new_set]))
    kernel_validation <- as.kernelMatrix(as.matrix(kernel_train[validation_set,train_new_set]))
    reg_ksvm <- ksvm(kernel_train_new,as.matrix(label_train),kernel='matrix',epsilon=parameter[i,2],type="nu-svr",cross=0,C = parameter[i,3],nu = parameter[i,4])
    sv_index <- reg_ksvm@alphaindex
    score_validation<- predict(reg_ksvm,as.kernelMatrix(kernel_validation[,sv_index]))
    cross_validation_error[1,j,2*i-1] <- reg_ksvm@error
    cross_validation_error[2,j,2*i-1] <-sum((score_validation-label_validation)^2)/length(label_validation)
    cross_validation_error[3,j,2*i-1] <- our_mean_error(reg_ksvm@fitted,ch1_train[train_new_set,c(1,14,12)])
    cross_validation_error[4,j,2*i-1] <- our_mean_error(score_validation,ch1_train[validation_set,c(1,14,12)])
    temp_train <- final_error(label_train[,1],as.numeric(reg_ksvm@fitted),ch1_train[train_new_set,c(1,14,12)])
    temp_test <-  final_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
    cross_validation_error[5,j,2*i-1] <- temp_train[1]
    cross_validation_error[6,j,2*i-1] <- temp_test[1]
    cross_validation_error[7,j,2*i-1] <- temp_train[2]
    cross_validation_error[8,j,2*i-1] <- temp_test[2]
    cross_validation_error[9,j,2*i-1] <- temp_train[3]
    cross_validation_error[10,j,2*i-1] <- temp_test[3]
    weak_leaner_matrix_train[,2*i-1] <- reg_ksvm@fitted
    weak_leaner_matrix_validation[,2*i-1] <- score_validation
    
    
    lambda <- parameter[i,5]
    ########
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
    ###############################
    count[] <- length(s)/14
    d <- 1/count[s[1:length(s)]]
    sqrt_D <- sqrt(length(s)) * diag(sqrt(d))/sqrt(sum(d))
    ##########
    score_validation <- kernel_validation %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
    score_validation_train <- kernel_train_new %*% sqrt_D %*% solve(lambda * diag(nrow(kernel_train_new)) + sqrt_D %*% kernel_train_new %*% sqrt_D) %*% sqrt_D %*% label_train
    cross_validation_error[1,j,2*i] <-sum((score_validation_train-label_train)^2)/length(label_train)
    cross_validation_error[2,j,2*i] <-sum((score_validation-label_validation)^2)/length(label_validation)
    cross_validation_error[3,j,2*i] <- our_mean_error(score_validation_train,ch1_train[train_new_set,c(1,14,12)])
    cross_validation_error[4,j,2*i] <- our_mean_error(score_validation,ch1_train[validation_set,c(1,14,12)])
    temp_train <- final_error(label_train[,1],as.numeric(score_validation_train),ch1_train[train_new_set,c(1,14,12)])
    temp_test <-  final_error(label_validation[,1],as.numeric(score_validation),ch1_train[validation_set,c(1,14,12)])
    cross_validation_error[5,j,2*i] <- temp_train[1]
    cross_validation_error[6,j,2*i] <- temp_test[1]
    cross_validation_error[7,j,2*i] <- temp_train[2]
    cross_validation_error[8,j,2*i] <- temp_test[2]
    cross_validation_error[9,j,2*i] <- temp_train[3]
    cross_validation_error[10,j,2*i] <- temp_test[3]
    weak_leaner_matrix_train[,2*i] <- score_validation_train
    weak_leaner_matrix_validation[,2*i] <- score_validation
    #}
  }
  
  ############################### glm
  ensemble_train <- as.data.frame(cbind(weak_leaner_matrix_train,label_train))
  f <- as.formula(paste("label_train ~", paste(name_col[!name_col %in% "label_train"], collapse = " + ")))
  glm_model <- glm(f, data = ensemble_train,family = gaussian(link = "identity"))
  train_predicted_glm <- predict(glm_model,as.data.frame(weak_leaner_matrix_train))
  cross_validation_error[1,j,2*i+1] <-sum(abs(train_predicted_glm-label_train)^2)/nrow(ensemble_train)
  cross_validation_error[9,j,2*i+1] <- temp_train[3]
  cross_validation_error[10,j,2*i+1] <- temp_test[3]
  ############################### PCR
  PCR_weight <- our_PCR(as.matrix(weak_leaner_matrix_train),label_train)
  PCR_score_train <- t(PCR_weight %*% t(weak_leaner_matrix_train))
  PCR_score_validation <- t(PCR_weight %*% t(weak_leaner_matrix_validation))
  cross_validation_error[1,j,2*i+2]<-sum(abs(PCR_score_train-label_train)^2)/nrow(ensemble_train)
  cross_validation_error[2,j,2*i+2] <-sum(abs(PCR_score_validation-label_validation)^2)/nrow(weak_leaner_matrix_validation)
  cross_validation_error[3,j,2*i+2] <- our_mean_error(PCR_score_train,ch1_train[train_new_set,c(1,14,12)])
  cross_validation_error[4,j,2*i+2] <- our_mean_error(PCR_score_validation,ch1_train[validation_set,c(1,14,12)])
  temp_train <- final_error(label_train[,1],as.numeric(PCR_score_train),ch1_train[train_new_set,c(1,14,12)])
  temp_test <-  final_error(label_validation[,1],as.numeric(PCR_score_validation),ch1_train[validation_set,c(1,14,12)])
  cross_validation_error[5,j,2*i+2] <- temp_train[1]
  cross_validation_error[6,j,2*i+2] <- temp_test[1]
  cross_validation_error[7,j,2*i+2] <- temp_train[2]
  cross_validation_error[8,j,2*i+2] <- temp_test[2]
  cross_validation_error[9,j,2*i+2] <- temp_train[3]
  cross_validation_error[10,j,2*i+2] <- temp_test[3]
  ################################# ridge
  ens_lamda <- 0.001
  ridge_score_validation <- weak_leaner_matrix_validation %*% solve(t(weak_leaner_matrix_train) %*% weak_leaner_matrix_train+ ens_lamda * diag(ncol(weak_leaner_matrix_train))) %*% t(weak_leaner_matrix_train) %*% label_train
  ridge_score_train <- weak_leaner_matrix_train %*% solve(t(weak_leaner_matrix_train) %*% weak_leaner_matrix_train+ ens_lamda * diag(ncol(weak_leaner_matrix_train))) %*% t(weak_leaner_matrix_train) %*% label_train
  cross_validation_error[1,j,2*i+3]<-sum(abs(ridge_score_train-label_train)^2)/nrow(ensemble_train)
  cross_validation_error[2,j,2*i+3] <-sum(abs(ridge_score_validation-label_validation)^2)/nrow(weak_leaner_matrix_validation)
  cross_validation_error[3,j,2*i+3] <- our_mean_error(ridge_score_train,ch1_train[train_new_set,c(1,14,12)])
  cross_validation_error[4,j,2*i+3] <- our_mean_error(ridge_score_validation,ch1_train[validation_set,c(1,14,12)])
  temp_train <- final_error(label_train[,1],as.numeric(ridge_score_train),ch1_train[train_new_set,c(1,14,12)])
  temp_test <-  final_error(label_validation[,1],as.numeric(ridge_score_validation),ch1_train[validation_set,c(1,14,12)])
  cross_validation_error[5,j,2*i+3] <- temp_train[1]
  cross_validation_error[6,j,2*i+3] <- temp_test[1]
  cross_validation_error[7,j,2*i+3] <- temp_train[2]
  cross_validation_error[8,j,2*i+3] <- temp_test[2]
  cross_validation_error[9,j,2*i+3] <- temp_train[3]
  cross_validation_error[10,j,2*i+3] <- temp_test[3]
  
}
#!/usr/bin/env Rscript


## libraries
library(kernlab)
library(plyr)
#library(deepnet)
# Arguments
#args=commandArgs(trailingOnly = TRUE)
############################################################# function
  #my_nn <- nn.train(x = cbind(weak_leaner_matrix_train,rep(1,nrow(weak_leaner_matrix_train))),y = label_train, hidden = c(20,5),activationfun = "sigm",
  #                  learningrate = 0.8, momentum = 0.4, learningrate_scale = 1, output = "linear",
  #                  numepochs = 1000, batchsize = 100, hidden_dropout = 0.3, visible_dropout = 0)
  #nn_score_train <- nn.predict(my_nn,cbind(weak_leaner_matrix_train,rep(1,nrow(weak_leaner_matrix_train))))
  #nn_score_validation <- nn.predict(my_nn,cbind(weak_leaner_matrix_validation,rep(1,nrow(weak_leaner_matrix_validation))))  
  #cross_validation_error[1,j,2*i+4]<-sum(abs(nn_score_train-label_train)^2)/nrow(ensemble_train)
  #cross_validation_error[2,j,2*i+4] <-sum(abs(nn_score_validation-label_validation)^2)/nrow(weak_leaner_matrix_validation)
  #cross_validation_error[3,j,2*i+4] <- our_mean_error(nn_score_train,ch1_train[train_new_set,c(1,14,12)])
  #cross_validation_error[4,j,2*i+4] <- our_mean_error(nn_score_validation,ch1_train[validation_set,c(1,14,12)])
  #cross_validation_error[5,j,2*i+4] <- global_error(label_train[,1],as.numeric(nn_score_train),ch1_train[train_new_set,c(1,14,12)])
  #cross_validation_error[6,j,2*i+4] <- global_error(label_validation[,1],as.numeric(nn_score_validation),ch1_train[validation_set,c(1,14,12)])
  
  
  ################################### GenSA
  #param_lower <- rep(0,ncol(weak_leaner_matrix_train))
  #param_upper <- rep(1,ncol(weak_leaner_matrix_train))
  #set.seed(1234)
  #out <- GenSA(lower = param_lower, upper = param_upper, fn = my_sa,
  #             control=list(maxit = 100,verbose=TRUE))
}
#write.table(cross_validation_error,"cv_error_sum.txt",quote = F,sep = "\t",row.names = F, col.names = F)
