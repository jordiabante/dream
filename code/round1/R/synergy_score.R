#!/usr/bin/env Rscript

## using library(kernlab
library(kernlab)

synergy_score <- function(kernel_train,kernel_test,label){
 
  
  ## random partition
  #random_set <- sample(1:nrow(kernel_train), nrow(kernel_train), replace = FALSE)
  #num_train = 2000;
  #train_index <-  random_set[1:num_train]
  #test_index <- random_set[2001:nrow(kernel_train)]
  
  ##filtering 
  #max_our =40
  #label[which(label > max_our )] <- max_our
  #label[which(label < -max_our )] <- -max_our
  #k_train <- as.kernelMatrix(as.matrix(kernel_train[train_index,train_index]))
  #k_test <- as.kernelMatrix(as.matrix(kernel_train[test_index,train_index]))
  #reg <- ksvm(k_train,as.matrix(label[train_index]),kernel='matrix',epsilon=1,type="nu-svr",cross=10,C = 50,nu = 1)
  
  
  ## ksvm
  k_train <- as.kernelMatrix(as.matrix(kernel_train))
  k_test <- as.kernelMatrix(as.matrix(kernel_test))
  reg <- ksvm(k_train,as.matrix(label),kernel='matrix',epsilon=1,type="nu-svr",cross=10,C = 50,nu = 1)
  sv_index <- reg@alphaindex
  score<- predict(reg,as.kernelMatrix(k_test[,sv_index]))
  
  #plot(label[test_index],score, xlim = c(-max_our,max_our), ylim = c(-max_our,max_our))

  write.csv(score,"score_gene_exp(all data).csv")
  ## 
  return(score)
  
}


kernel_function <- function(x,y){
  index_cell_x <- which(x[1] == cell_line_order)
  index_drug_comb_x <- which(x[14] == drug_combination_order)
  index_cell_y <- which(y[1] == cell_line_order)
  index_drug_comb_y <- which(y[14] == drug_combination_order)
  return (kernel_cell_line[index_cell_x,index_cell_y] +  kernel_drug[index_drug_comb_x,index_drug_comb_y])
}

