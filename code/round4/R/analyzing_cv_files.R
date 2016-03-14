Dir_path_svr <- "../../data/round4/cv/svr"
fi_svr<-list.files(Dir_path_svr,full.names=T)
dat_svr<-lapply(fi_svr,read.table,sep="\t",header=F)
num_week_learner <- length(fi_svr)
svr_num_param <- 3
Dir_path_ridge <- "../../data/round4/cv/ridge"
fi_ridge<-list.files(Dir_path_ridge,full.names=T)
dat_ridge<-lapply(fi_ridge,read.table,sep="\t",header=F)
num_week_learner <- length(fi_ridge)
ridge_num_param <- 1
param_global_product <- matrix(0,nrow = num_week_learner,ncol = svr_num_param+ridge_num_param+1)
for(i in 1:num_week_learner){
  temp_matrix <- as.matrix(dat_svr[[i]])
  param_global_product[i,1] <- substr(fi_svr[i],nchar(Dir_path_svr)+9,nchar(fi_svr[i]))
  param_global_product[i,2:4] <- temp_matrix[which.max(temp_matrix[,4]),1:3]
  temp_matrix <- as.matrix(dat_ridge[[i]])
  param_global_product[i,5] <- temp_matrix[which.max(temp_matrix[,2]),1]
}
write.table(param_global_product,"../../data/round4/param/param_global_product_non_weighted.txt",sep = "\t",quote = F,row.names = F,col.names = F)
