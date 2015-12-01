#!/usr/bin/env Rscript

# Libraries
library(Matrix)

# Arguments
args=commandArgs(trailingOnly = TRUE)
leaderboard_test_path=args[1]
leaderboard_train_path=args[2]
kernel_cell_line_path=args[3]
kernel_drug_path=args[4]
cell_line_order_path=args[5]
drug_combination_order_path=args[6]
outsuffix=args[7]

# Read data
leaderboard_test=read.table(leaderboard_test_path,sep=",",header=T)
leaderboard_train=read.table(leaderboard_train_path,sep=",",header=T)
kernel_cell_line=read.table(kernel_cell_line_path,sep=",",header=T)
kernel_drug=read.table(kernel_drug_path,sep=",",header=T)
cell_line_order=read.table(cell_line_order_path,sep=",")
drug_combination_order=read.table(drug_combination_order_path,sep=",")

# Functions
make_total_kernel <- function(leaderboard_table_test,leaderboard_table_train){

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
  kernel_total_train <- as.matrix(forceSymmetric(as.matrix(kernel_total)))

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

  # Write output
  write.table(kernel_total_train,file=paste("../../data/round1/kernels/kernel_train_",outsuffix,".csv"),col.names=F,row.names=F,sep=",",quote=F)
  write.table(kernel_test,file=paste("../../data/round1/kernels/kernel_test_",outsuffix,".csv"),col.names=F,row.names=F,sep="\t",quote=F)
}

# Call function
make_total_kernel(leaderboard_test,leaderboard_train)
