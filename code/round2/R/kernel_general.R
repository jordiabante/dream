#!/usr/bin/env Rscript

# Libraries
library(Matrix)

# Arguments
#args=commandArgs(trailingOnly = TRUE)
#leaderboard_test_path=args[1]
#leaderboard_train_path=args[2]
#kernel_cell_line_path=args[3]
#kernel_drug_path=args[4]
#cell_line_order_path=args[5]
#drug_combination_order_path=args[6]
#outsuffix=args[7]

# Read data
#leaderboard_test=read.table(leaderboard_test_path,sep=",",header=T)
#leaderboard_train=read.table(leaderboard_train_path,sep=",",header=T)
#kernel_cell_line=read.table(kernel_cell_line_path,sep="\t",header=F)
#kernel_drug=read.table(kernel_drug_path,sep=",",header=T)
#cell_line_order=read.table(cell_line_order_path,sep=",")
#drug_combination_order=read.table(drug_combination_order_path,sep=",")

leaderboard_test=read.table("ch1_leaderBoard_monoTherapy.csv",sep=",",header=T)
leaderboard_train=read.table("ch1_train_combination_and_monoTherapy.csv",sep=",",header=T)
kernel_cell_line=read.table("dot_product_mutations.csv",sep=",",header=F)
kernel_drug=read.table("drug_kernel_dot.csv",sep=",",header=F)
cell_line_order=as.data.frame(read.table("cell_line_order.csv",sep=","),header=F)
drug_combination_order=read.table("drug_comb_name.csv",sep=",",header=F)



header_information <- names(leaderboard_test)
drug_index <- which(header_information == "COMBINATION_ID")
cell_index <- which(header_information == "CELL_LINE")
quality_index <- which(header_information == "QA")
synergy_index <- which(header_information == "SYNERGY_SCORE")
leaderboard_train <- leaderboard_train[which(leaderboard_train[,quality_index]==1),]
train_score <- as.data.frame(leaderboard_train[,synergy_index])

############################################################ should change ###############
write.table(train_score,"train_score.csv",sep=",",col.names=F,row.names=F,quote=F)
#########################################################################################


kernel_total <- matrix(0,nrow = nrow(leaderboard_train),ncol = nrow(leaderboard_train))
for (i in 1:nrow(kernel_total)){
  for (j in i:nrow(kernel_total)){
    kernel_total[i,j] <- (kernel_cell_line[as.integer(leaderboard_train[i,cell_index]),as.integer(leaderboard_train[j,cell_index])] +  kernel_drug[as.integer(leaderboard_train[i,drug_index]),as.integer(leaderboard_train[j,drug_index])])/2
  }
}
kernel_total_train <- as.matrix(forceSymmetric(as.matrix(kernel_total)))

## finding kernel matrix for test data
kernel_test <- matrix(0,nrow = nrow(leaderboard_test),ncol = nrow(leaderboard_train))
for(i in 1:nrow(kernel_test)){
  for (j in 1:ncol(kernel_test)){
    kernel_test[i,j] <- (kernel_cell_line[as.integer(leaderboard_test[i,cell_index]),as.integer(leaderboard_train[j,cell_index])] +  kernel_drug[as.integer(leaderboard_test[i,drug_index]),as.integer(leaderboard_train[j,drug_index])])/2
  }
}

# Write output
#write.table(kernel_total_train,file=paste("../../data/round1/kernels/kernel_train_",outsuffix,".csv"),col.names=F,row.names=F,sep=",",quote=F)
#write.table(kernel_test,file=paste("../../data/round1/kernels/kernel_test_",outsuffix,".csv"),col.names=F,row.names=F,sep="\t",quote=F)

write.table(kernel_total_train,file="kernel_train_mutation.csv",col.names=F,row.names=F,sep=",",quote=F)
write.table(kernel_test,file="kernel_test_mutation.csv",col.names=F,row.names=F,sep="\t",quote=F)
