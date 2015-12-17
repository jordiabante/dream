#!/usr/bin/env Rscript

# Libraries
library(Matrix)

########################################################## function ##################################################
make_kernel <- function(kernel_cell_line,kernel_drug,outsuffix){  
  header_information <- names(leaderboard_test)
  drug_index <- which(header_information == "COMBINATION_ID")
  cell_index <- which(header_information == "CELL_LINE")
  quality_index <- which(header_information == "QA")
  synergy_index <- which(header_information == "SYNERGY_SCORE")
  leaderboard_train <- leaderboard_train[which(leaderboard_train[,quality_index]==1),]
  train_score <- as.data.frame(leaderboard_train[,synergy_index])
  
  # write.table(train_score,"train_score.csv",sep=",",col.names=F,row.names=F,quote=F)
  # write.table(train_score,"../../data/originals/train_score.csv",sep=",",col.names=F,row.names=F,quote=F)
  
  
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
  write.table(kernel_total_train,file=paste("../../data/round2/kernel_train_test/kernel_train_",outsuffix,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  write.table(kernel_test,file=paste("../../data/round2/kernel_train_test/kernel_test_",outsuffix,".txt",sep=""),col.names=F,row.names=F,sep="\t",quote=F)
  
  # write.table(kernel_total_train,file="kernel_train_mutation.csv",col.names=F,row.names=F,sep=",",quote=F)
  #write.table(kernel_test,file="kernel_test_mutation.csv",col.names=F,row.names=F,sep="\t",quote=F)
}

##################################################################################################################################################################
# Read data
leaderboard_test=read.table("../../data/originals/ch1_leaderBoard_monoTherapy.csv",sep=",",header=T)
leaderboard_train=read.table("../../data/originals/ch1_train_combination_and_monoTherapy.csv",sep=",",header=T)
cell_line_order=read.table("../../data/originals/cell_line_order.csv",sep=",",,header=F)
drug_combination_order=read.table("../../data/originals/drug_comb_name.csv",sep=",",,header=F)
kernel_names=read.table("../../data/round2/kernels/kernel_names.txt",sep="\t",header=T)

#leaderboard_test=read.table("ch1_leaderBoard_monoTherapy.csv",sep=",",header=T)
#leaderboard_train=read.table("ch1_train_combination_and_monoTherapy.csv",sep=",",header=T)
#cell_line_order=as.data.frame(read.table("cell_line_order.csv",sep=","),header=F)
#drug_combination_order=read.table("drug_comb_name.csv",sep=",",header=F)
#kernel_names=read.table("kernel_names.txt",sep="\t",header=T)

name <- names(kernel_names)
drug_index <- which((name=="drug_dot_product") | (name=="drug_angular_similarity")  )
cell_index <- which(!(name=="drug_dot_product") & !(name=="drug_angular_similarity"))

for(i in  1:length(drug_index)){
  for(j in 1:nrow(kernel_names)){
    kernel_drug_path <- kernel_names[j,drug_index[i]]
    if(!is.na(kernel_drug_path)){
      drug_outsuffix <- sub('../../data/round2/kernels/','',kernel_drug_path)
      kernel_drug <- read.table(paste(kernel_drug_path,".txt",sep=""),sep="\t",header=F)      
      for(l in 1:length(cell_index)){
        if((l==3) | (l==4) | (l==7) | (l==8)){
          for(k in 1:nrow(kernel_names)){
            kernel_cell_path <- kernel_names[k,cell_index[l]]
            if(!is.na(kernel_cell_path)){
              cell_outsuffix <-sub('../../data/round2/kernels/','',kernel_cell_path)
              kernel_cell_line <- read.table(paste(kernel_cell_path,".txt",sep=""),sep="\t",header=F)
              print(paste("kernel train and test : ",drug_outsuffix,"_",cell_outsuffix," has been started",sep=""))
              make_kernel(kernel_cell_line,kernel_drug,paste(drug_outsuffix,"_",cell_outsuffix,sep=""))
              print(paste("kernel train and test : ",drug_outsuffix,"_",cell_outsuffix," has been completed",sep=""))
            }	  
          }
        }
      }
    }
  }
}





