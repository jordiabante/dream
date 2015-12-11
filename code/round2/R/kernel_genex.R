#!/usr/bin/env Rscript
################################
## Data set:
## 83 cell lines, 17419 genes
## Missing cell lines: MDA-MB-175-VII and NCI-H1437
################################

## Dependancies
#require('SpatialEpi')
#require('kernlab')

# Read raw data
read.table('../../data/originals/gex.csv.gz',sep=",",header=T,row.names=1)->x
############################################# Correlation ##########################################
## No need to normalize
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
write.table(cls_corr,file="../../data/round2/kernels/corr_genex.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

########################################### Kernel b.f. ########################################
####################### Dot product before filtering
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

write.table(dot_product,file="../../data/round2/kernels/dot_product_genex_original.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

####################### Angular similarity
angular_similarity=matrix(0,nrow=ncol(dot_product),ncol=ncol(dot_product))
for(i in 1:nrow(dot_product))
{
    for(j in 1:nrow(dot_product))
    {
        if(i==j)
        {
            angular_similarity[i,j]=1
        } else {
            angular_similarity[i,j]=1-acos(dot_product[i,j])/3.14159
        }
    }
}
write.table(angular_similarity,file="../../data/round2/kernels/angular_similarity_genex_original.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

############################################# Kernel a.f. #############################################
####################### Dot product after filtering
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
indexes=variances[1:(0.01*nrow(variances)),2]
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
write.table(dot_product,file="../../data/round2/kernels/dot_product_genex_filtered.txt",
            col.names=F,row.names=F,sep="\t",quote=F)

####################### Angular similarity
angular_similarity=matrix(0,nrow=ncol(dot_product),ncol=ncol(dot_product))
for(i in 1:nrow(dot_product))
{
    for(j in 1:nrow(dot_product))
    {
        if(i==j)
        {
            angular_similarity[i,j]=1
        } else {
            angular_similarity[i,j]=1-acos(dot_product[i,j])/3.14159
        }
    }
}
write.table(angular_similarity,file="../../data/round2/kernels/angular_similarity_genex_filtered.txt",
            col.names=F,row.names=F,sep="\t",quote=F)
print("kernel genex completed")
