#!/usr/bin/env Rscript


############################################# kernel functions ###################################

My_rbf <- function(cell_level,feature,sigma){
  rbf_kernel <- matrix(0,ncol = nrow(cell_level), nrow = nrow(cell_level))
  for (i in 1:nrow(cell_level)) {
    for (j in 1:nrow(cell_level)) { 
      rbf_kernel[i,j] <- exp(-sigma * norm(feature[i,]-feature[j,],"2")^2)
    }
  }
  return(rbf_kernel)
}


inner_angular_product <- function(cell_level,feature){
  my_kernel_inner <- matrix(0,nrow = nrow(cell_level),ncol = nrow(cell_level))
  my_kernel_angular <- matrix(0,nrow = nrow(cell_level),ncol = nrow(cell_level))
  for (i in 1:nrow(cell_level)){
    for (j in 1:nrow(cell_level)){
      my_kernel_inner[i,j] <- (feature[i,] %*% feature[j,])/(norm(feature[i,],"2")*norm(feature[j,],"2"))
      my_kernel_angular[i,j] <- 1-2*acos(pmin(pmax(my_kernel_inner[i,j],-1.0),1.0))/pi
    }
  }
  return(rbind(my_kernel_inner,my_kernel_angular))
}



 
######################################## Reading file and Initialization ###########################
#print("Start reading Files and Initializations")
mutation=read.table('../../data/originals/mutations.csv.gz',sep=",",header=T)
#mutation <- read.table('mutations.csv.gz',sep=",",header=T)
 header_information <- names(mutation)
 gene_index <- which(header_information == "Gene.name")
 cell_index <- which(header_information == "cell_line_name")
 cell_level <- as.data.frame(levels(mutation[,cell_index]))
 Fathmm_index <- which(header_information=="FATHMM.prediction")
 description_index <- which(header_information=="Mutation.Description")
 
 # romving "_" in data
 mutation[,gene_index] <- as.data.frame.character(sapply(mutation[,gene_index], function(x) strsplit(as.character(x), split='_', fixed=TRUE)[[1]][1]))
 
 # extracting all genes
 gene_level <- as.data.frame(levels(mutation[,gene_index]))
 
 # extracting driver genes
 cancer_gene <- as.data.frame(mutation[which(mutation[,Fathmm_index]=="CANCER"),gene_index])
 cancer_gene <- as.data.frame(cancer_gene[!duplicated(cancer_gene),1 ])
 
 # extracting non-silent mutations
 mutation_non_silent <- mutation[which(mutation[,description_index]!="Substitution - coding silent"),]
 gene_level_nonsilent <- as.data.frame(levels(mutation_non_silent[,gene_index]))

 
############################################# constructing features #################################### 
#print("start constructing features") 
# feature all genes
 feature_all_genes <- matrix(0,nrow = nrow(cell_level),ncol = nrow(gene_level))
  for (i in 1:nrow(cell_level)){
      index <- which(mutation[,cell_index]==cell_level[i,1])
      for(j in 1:length(index)){
        feature_all_genes[i,as.integer(mutation[index[j],gene_index])] <- 1
      }
  }

 # feature driver genes
  feature_driver_genes <- matrix(0,nrow = nrow(cell_level),ncol = nrow(cancer_gene))
  for(i in 1:nrow(cancer_gene)){
    index <- which(mutation[,gene_index]==cancer_gene[i,1])
    if(length(index)>0){
      for(j in 1:length(index)){
        feature_driver_genes[as.integer(mutation[index[j],cell_index]),i] <- 1
      }
    }
  }

# nonsilent driver mutations
  feature_nonsilent <- matrix(0,nrow = nrow(cell_level),ncol = nrow(gene_level_nonsilent))
  for (i in 1:nrow(cell_level)){
    index <- which(mutation_non_silent[,cell_index]==cell_level[i,1])
    for(j in 1:length(index)){
      feature_nonsilent[i,as.integer(mutation_non_silent[index[j],gene_index])] <- 1
    }
  }
  
################################################### constructing kernels ##########################  
#print("start constructing kernels")  
  
  # Kernel inner and angular product
  temp <- inner_angular_product(cell_level,feature_all_genes)
  kernel_all_genes_dot <- temp[1:nrow(cell_level),]
  kernel_all_genes_angular <- temp[(nrow(cell_level)+1):(2*nrow(cell_level)),]
  
  temp <- inner_angular_product(cell_level,feature_driver_genes)
  kernel_driver_genes_dot <- temp[1:nrow(cell_level),]
  kernel_driver_genes_angular <- temp[(nrow(cell_level)+1):(2*nrow(cell_level)),]
  
  temp <- inner_angular_product(cell_level,feature_nonsilent)
  kernel_nonsilent_dot <- temp[1:nrow(cell_level),]
  kernel_nonsilent_angular <- temp[(nrow(cell_level)+1):(2*nrow(cell_level)),]
  
  
  # kernel rbf
  sigma <- 0.005
  kernel_all_genes_rbf <- My_rbf(cell_level,feature_all_genes,sigma)
  kernel_driver_genes_rbf <- My_rbf(cell_level,feature_driver_genes,sigma)
  kernel_nonsilent_rbf <- My_rbf(cell_level,feature_nonsilent,sigma)
  
  # kernel correlation
  kernel_all_genes_corr <- cor(t(feature_all_genes))
  kernel_driver_genes_corr <- cor(t(feature_driver_genes))
  kernel_nonsilent_corr <- cor(t(feature_nonsilent))
  
################################################# saving files ###################################  
  write.table(kernel_all_genes_dot,file="../../data/round2/kernels/dot_product_mutations_original.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(kernel_all_genes_angular,file="../../data/round2/kernels/angular_similarity_mutations_original.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 
  
  write.table(kernel_driver_genes_dot,file="../../data/round2/kernels/dot_product_mutations_driver.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(kernel_driver_genes_angular,file="../../data/round2/kernels/angular_similarity_mutations_driver.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 
  write.table(kernel_nonsilent_dot,file="../../data/round2/kernels/dot_product_mutations_nonsilent.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(kernel_nonsilent_angular,file="../../data/round2/kernels/angular_similarity_mutations_nonsilent.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 
  
  write.table(kernel_all_genes_corr,file="../../data/round2/kernels/corr_mutations_original.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(kernel_driver_genes_corr,file="../../data/round2/kernels/corr_mutations_driver.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 
  write.table(kernel_nonsilent_corr,file="../../data/round2/kernels/corr_mutations_nonsilent.txt",
              sep="\t",quote=F,col.names=F,row.names=F)

  write.table(kernel_all_genes_rbf,file="../../data/round2/kernels/rbf_mutations_original.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
  write.table(kernel_driver_genes_rbf,file="../../data/round2/kernels/rbf_mutations_driver.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 
  write.table(kernel_nonsilent_rbf,file="../../data/round2/kernels/rbf_mutations_nonsilent.txt",
              sep="\t",quote=F,col.names=F,row.names=F)

 write.table(feature_all_genes,file="../../data/round2/features/mutation_features_original.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
 write.table(feature_driver_genes,file="../../data/round2/features/mutation_features_driver.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
  
 write.table(feature_nonsilent,file="../../data/round2/features/mutation_features_nonsilent.txt",
              sep="\t",quote=F,col.names=F,row.names=F)
print("mutation kernels completed")
