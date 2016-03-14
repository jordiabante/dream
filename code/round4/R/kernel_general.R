#!/usr/bin/env Rscript

# Libraries
library(Matrix)

########################################################## function ##################################################
make_kernel <- function(kernel_cell_line,kernel_drug,outsuffix){  
  header_information <- names(test_data)
  drug_index <- which(header_information == "COMBINATION_ID")
  cell_index <- which(header_information == "CELL_LINE")
  quality_index <- which(header_information == "QA")
  synergy_index <- which(header_information == "SYNERGY_SCORE")
  leaderboard_train <- leaderboard_train[which(leaderboard_train[,quality_index]==1),]
  
  
  kernel_total <- matrix(0,nrow = nrow(leaderboard_train),ncol = nrow(leaderboard_train))
  for (i in 1:nrow(kernel_total)){
    for (j in i:nrow(kernel_total)){
      kernel_total[i,j] <- (kernel_cell_line[as.integer(leaderboard_train[i,cell_index]),as.integer(leaderboard_train[j,cell_index])] *  kernel_drug[as.integer(leaderboard_train[i,drug_index]),as.integer(leaderboard_train[j,drug_index])])
    }
  }
  kernel_total_train <- as.matrix(forceSymmetric(as.matrix(kernel_total)))
  
  ## finding kernel matrix for test data
  kernel_test <- matrix(0,nrow = nrow(leaderboard_test),ncol = nrow(leaderboard_train))
  for(i in 1:nrow(kernel_test)){
    for (j in 1:ncol(kernel_test)){
      kernel_test[i,j] <- (kernel_cell_line[as.integer(leaderboard_test[i,cell_index]),as.integer(leaderboard_train[j,cell_index])] *  kernel_drug[as.integer(leaderboard_test[i,drug_index]),as.integer(leaderboard_train[j,drug_index])])
    }
  }
  
  ## finding kernel matrix for test data
  kernel_test_final <- matrix(0,nrow = nrow(test_data),ncol = nrow(leaderboard_train))
  for(i in 1:nrow(kernel_test_final)){
    for (j in 1:ncol(kernel_test_final)){
      kernel_test_final[i,j] <- (kernel_cell_line[as.integer(test_data[i,cell_index]),as.integer(leaderboard_train[j,cell_index])] *  kernel_drug[as.integer(test_data[i,drug_index]),as.integer(leaderboard_train[j,drug_index])])
    }
  }
  # Write output
  write.table(kernel_total_train,file=paste("../../data/round4/kernel_train_test_product/kernel_train_",outsuffix,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(kernel_test,file=paste("../../data/round4/kernel_train_test_product/kernel_test_",outsuffix,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(kernel_test_final,file=paste("../../data/round4/kernel_train_test_product/kernel_final_test_",outsuffix,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
}


# Read data
leaderboard_test=read.table("../../data/originals/ch1_leaderBoard_monoTherapy.csv",sep=",",header=T)
leaderboard_train=read.table("../../data/originals/ch1_train_combination_and_monoTherapy.csv",sep=",",header=T)
cell_line_order=read.table("../../data/originals/cell_line_order.csv",sep=",",,header=F)
drug_combination_order=read.table("../../data/originals/drug_comb_order.csv",sep=",",,header=F)

Dir_path <- "../../data/round4/kernels"
fi<-list.files(Dir_path,full.names=T)
fi <- fi[-which(grepl("rbf",fi)==TRUE)]
num_matrix <- length(fi)
drug_index <- which(grepl("drug",fi)==TRUE)
cell_index <- which(!(grepl("drug",fi)==TRUE))

for(i in drug_index){
  kernel_drug_path <- fi[i]
  drug_outsuffix <- sub(paste(Dir_path,"/",sep = ""),'',kernel_drug_path)
  drug_outsuffix <-sub(".txt",'',drug_outsuffix)
  kernel_drug <- read.table(kernel_drug_path,sep="\t",header=F)  
  print(drug_outsuffix)
  for(l in cell_index){
    kernel_cell_path <- fi[l]
    cell_outsuffix <-sub(paste(Dir_path,"/",sep = ""),'',kernel_cell_path)
    cell_outsuffix <-sub(".txt",'',cell_outsuffix)
    kernel_cell_line <- read.table(kernel_cell_path,sep="\t",header=F)
    print(cell_outsuffix)
    make_kernel(kernel_cell_line,kernel_drug,paste(drug_outsuffix,"_",cell_outsuffix,sep=""))
  }
}




