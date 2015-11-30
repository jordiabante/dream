#!/usr/bin/env Rscript
################################
## Data set:
## 83 cell lines, 17419 genes
## Missing cell lines: MDA-MB-175-VII and NCI-H1437
################################

## Dependancies
require('SpatialEpi')
require('kernlab')

## Plain correlation
## No need to normalize
read.table('../../data/originals/gex.csv.gz',sep=",",header=T,row.names=1)->x
cor(x,method="spearman")->cls_corr

# Add MDA-MB-175-VII
aux=rbind(cls_corr[1:43,],rep(-1,83),cls_corr[44:83,])
aux=cbind(aux[,1:43],rep(-1,84),aux[,44:83])
aux[44,44]=1
colnames(aux)[44]="MDA-MB-175-VII"
rownames(aux)[44]="MDA-MB-175-VII"
cls_corr=aux
# Add NCI-H1437
aux=rbind(cls_corr[1:52,],rep(-1,84),cls_corr[53:84,])
aux=cbind(aux[,1:52],rep(-1,85),aux[,53:84])
aux[53,53]=1
colnames(aux)[53]="NCI-H1437"
rownames(aux)[53]="NCI-H1437"
cls_corr=aux
write.table(cls_corr,file="../../data/round1/kernels/corr_genes.tsv",
            col.names=F,row.names=F,sep="\t",quote=F)

## Distance matrix
## Range is not between 0-1
as.matrix(dist(t(x),diag=T,upper=T,method="euclidean"))->cls_dist

## Kernel
## Dot product
dot_product=matrix(0,nrow=ncol(x),ncol=ncol(x))
colnames(dot_product)=colnames(x)
rownames(dot_product)=colnames(x)
for(i in 1:ncol(x))
{
    for(j in 1:ncol(x))
    {
        vector_product=x[,i]%*%x[,j]
        norm_first=sqrt(sum(x[,i]^2))
        norm_second=sqrt(sum(x[,j]^2))
        dot_product[i,j]=vector_product/(norm_first*norm_second)
    }
}
# Add MDA-MB-175-VII
aux=rbind(dot_product[1:43,],rep(-1,83),dot_product[44:83,])
aux=cbind(aux[,1:43],rep(-1,84),aux[,44:83])
aux[44,44]=1
colnames(aux)[44]="MDA-MB-175-VII"
rownames(aux)[44]="MDA-MB-175-VII"
dot_product=aux
# Add NCI-H1437
aux=rbind(dot_product[1:52,],rep(-1,84),dot_product[53:84,])
aux=cbind(aux[,1:52],rep(-1,85),aux[,53:84])
aux[53,53]=1
colnames(aux)[53]="NCI-H1437"
rownames(aux)[53]="NCI-H1437"
dot_product=aux
write.table(dot_product,file="../../data/round1/kernels/dot_product_genex.tsv",
            col.names=F,row.names=F,sep="\t",quote=F)

## Gaussian kernel
# Normalization
#x_norm=matrix(0,nrow=nrow(x),ncol=ncol(x))
#colnames(x_norm)=colnames(x)
#rownames(x_norm)=rownames(x)
#for(i in 1:ncol(x))
#{
#    x_norm[,i]=normalize(x[,i])
#}
## kernel
#gaussian_kernel=matrix(0,nrow=ncol(x),ncol=ncol(x))
#rbf=rbfdot(sigma = 1)
#colnames(gaussian_kernel)=colnames(x)
#rownames(gaussian_kernel)=colnames(x)
### Gaussian kernel sigma=1
#for(i in 1:nrow(gaussian_kernel))
#{
#    for(j in 1:ncol(gaussian_kernel))
#    {
#        gaussian_kernel[i,j]=rbf(x_norm[,i],x_norm[,j])
#    }   
#}


