#!/usr/bin/env Rscript


kernel_mutation <- function(mytable){
  header_information <- names(mytable)
  gene_name_index <- which(header_information == "Gene.name")
  cell_line_index <- which(header_information == "cell_line_name")
  
  my_data_frame <- data.frame(gene_names = mytable[,gene_name_index],cell_line_name = mytable[,cell_line_index],stringsAsFactors = TRUE)
  genes_level <- levels(my_data_frame[,1])
  cell_line_level <- levels(my_data_frame[,2])
  
  my_feature <- matrix(0,nrow = length(cell_line_level),ncol = length(genes_level))
  for (i in 1:length(cell_line_level)){
    index <- which(my_data_frame[,2] == cell_line_level[i])
    for (j in 1:length(index)){
      my_feature[i,as.integer(my_data_frame[index[j],1])] <- 1
    }
  }
  write.table(my_feature,file="../../data/round1/features/mutations_features.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  # inner product
  my_kernel <- matrix(0,nrow = length(cell_line_level),ncol = length(cell_line_level))
  for (i in 1:length(cell_line_level)){
    for (j in 1:length(cell_line_level)){
      my_kernel[i,j] <- (my_feature[i,] %*% my_feature[j,])/(norm(my_feature[i,],"2")*norm(my_feature[j,],"2"))
    }
  }
  write.table(my_kernel,file="../../data/round1/kernels/dot_product_mutations.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  # correlation
  my_kernel_2 <- cor(t(my_feature))
  
  write.table(my_kernel_2,file="../../data/round1/kernels/corr_mutations.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
}


mutation_data=read.table('../../data/originals/mutations.csv.gz',sep=",",header=T)
kernel_mutation(mutation_data)
