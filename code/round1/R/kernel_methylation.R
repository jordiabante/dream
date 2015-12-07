#!/usr/bin/env Rscript
####################################
## Data set:
## 84 cell lines, 104426 CpG shores
## Missing cell lines: SW620 (74)
####################################

## Dependancies
require('kernlab')

## Plain correlation
## No need to normalize
read.table('../../data/originals/methylation.csv.gz',sep=",",header=T,row.names=1,fill=TRUE)->x
x=x[complete.cases(x),]
cor(x,method="spearman")->cls_corr

# Add SW620
aux=rbind(cls_corr[1:73,],rep(-1,84),cls_corr[74:84,])
aux=cbind(aux[,1:73],rep(-1,85),aux[,74:84])
aux[74,74]=1
colnames(aux)[74]="SW620"
rownames(aux)[74]="SW620"
cls_corr=aux
write.table(cls_corr,file="../../data/round1/kernels/corr_methylation.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

## Kernel
# Get variance across cell lines for each gene
variances=matrix(0,nrow=nrow(x),ncol=2)
for(i in 1:nrow(x))
{
    variances[i,1]=var(as.numeric(x[i,]))
    variances[i,2]=rownames(x[i,])
}
# Rank genes based on variance
variances=variances[order(as.numeric(variances[,1]),decreasing=TRUE),]

# Get 10% genes with more differences
indexes=variances[1:(0.1*nrow(variances)),2]
x=x[indexes,]

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
# Add SW620
aux=rbind(dot_product[1:73,],rep(-1,84),dot_product[74:84,])
aux=cbind(aux[,1:73],rep(-1,85),aux[,74:84])
aux[74,74]=1
colnames(aux)[74]="SW620"
rownames(aux)[74]="SW620"
dot_product=aux
write.table(dot_product,file="../../data/round1/kernels/dot_product_methylation.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

