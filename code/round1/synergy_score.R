
synergy_score <- function(kernel_train,kernel_test,label){
 
  ## using library(kernlab
  library(kernlab)
  
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











make_total_kernel <- function(leaderboard_table_test,leaderboard_table_train){
  
  ## imbalanced issue
  #label <- leaderboard_table_train[,12]
  #pos_index <- which(label >0)
  #neg_index <- which(label <0)
  #neg_index_new <- sample(neg_index , length(pos_index)-length(neg_index),replace = TRUE) 
  
  #leaderboard_table_train<- rbind(leaderboard_table_train,leaderboard_table_train[neg_index_new,])
  
  #write(neg_index_new,"index_new_over.csv")
  
  
  kernel_total <- matrix(0,nrow = nrow(leaderboard_table_train),ncol = nrow(leaderboard_table_train))
  for (i in 1:nrow(kernel_total)){
    if (i %% 100 == 0){
      print(i)
    }
    index_cell_i <- which(leaderboard_table_train[i,1] == cell_line_order)
    index_drug_comb_i <- which(leaderboard_table_train[i,14] == drug_combination_order)
    for (j in i:ncol(kernel_total)){
      index_cell_j <- which(leaderboard_table_train[j,1] == cell_line_order)
      index_drug_comb_j <- which(leaderboard_table_train[j,14] == drug_combination_order)
      kernel_total[i,j] <- kernel_cell_line[index_cell_i,index_cell_j] +  kernel_drug[index_drug_comb_i,index_drug_comb_j]
    }
  }
  library(Matrix)
  kernel_total_train <- as.matrix(forceSymmetric(as.matrix(kernel_total))) 
  
  write.csv(kernel_total_train,"kernel_train_gene_exp(all data).csv")
  
  ## finding kernel matrix for test data
  kernel_test <- matrix(0,nrow = nrow(leaderboard_table_test),ncol = nrow(leaderboard_table_train))
  for(i in 1:nrow(kernel_test)){
    if (i %% 100 == 0){
      print(i)
    }
    index_cell_i <- which(leaderboard_table_test[i,1] == cell_line_order)
    temp <- as.character(leaderboard_table_test[i,14])
    index_drug_comb_i <- which(drug_combination_order == temp )
    for (j in 1:ncol(kernel_test)){
      index_cell_j <- which(leaderboard_table_train[j,1] == cell_line_order)
      index_drug_comb_j <- which(leaderboard_table_train[j,14] == drug_combination_order)
      kernel_test[i,j] <- kernel_cell_line[index_cell_i,index_cell_j] +  kernel_drug[index_drug_comb_i,index_drug_comb_j]
    }
  }
  write.csv(kernel_test,"kernel_test_gene_exp(all data).csv")
}


making_submission <- function(score,leaderboard_table_test){
  submission <-leaderboard_table_test[,c(1,14)]
  submission [,3] <-  score
  write.csv(submission,"submission_1B_filter.csv")
}
