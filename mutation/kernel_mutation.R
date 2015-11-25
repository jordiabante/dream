
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
  write.csv(my_feature,"mutation_feature (simple v.1).csv")
  
  # inner product
  my_kernel <- matrix(0,nrow = length(cell_line_level),ncol = length(cell_line_level))
  for (i in 1:length(cell_line_level)){
    for (j in 1:length(cell_line_level)){
      my_kernel[i,j] <- (my_feature[i,] %*% my_feature[j,])/(norm(my_feature[i,],"2")*norm(my_feature[j,],"2"))
    }
  }
  write.csv(my_kernel,"kernel_mutation (inner product)(simple v.1).csv")
  
  # correlation
  my_kernel_2 <- cor(t(my_feature))
  
  write.csv(my_kernel_2,"kernel_mutation (corr product)(simple v.1).csv")
}





