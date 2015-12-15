#!/usr/bin/env Rscript

library(org.Hs.eg.db)
library(KEGGREST)
library(stringr)
library(SpatialEpi)
library(kernlab)

## get the drug and target pairs
drug_info = read.csv(file='../../data/originals/drug_info_release.csv')
target = drug_info[,2]
drug = drug_info[,1]
# Delete lines that have "DNA" and "Methylation"
drug = drug[-grep("DNA", target)]
target = target[-grep("DNA", target)]
drug = as.character(drug[-grep("Methylation", target)])
target = as.character(target[-grep("Methylation", target)])

# extract gene information from target, as a class of list
genelist = rep(NA, length(target))
for (i in 1:length(target))
{
  target[i] = gsub("[ ]", "", target[i])
  genelist[i] = strsplit(target[i], split = ",")
}

## get the target and pathway pair
path_genelist = rep(NA, length(genelist))
for (i in 1:length(genelist))
{
  cat(i)
  # get ensemble ID for each gene list
  ID = paste(unlist(mget(unlist(genelist[i]), org.Hs.egSYMBOL2EG, ifnotfound = NA)))
  ID = paste("hsa:", ID)
  ID = gsub("[ ]", "", ID)
  pathway = paste(unlist(keggLink("pathway", ID)))
  
  # store the pathway into the path_genelist vector
  if(length(pathway) !=0)
  {
    pathway = sapply(strsplit(pathway, split=':', fixed=TRUE), function(x) (x[2]))
    pathway = unique(pathway)
    path_genelist[i] = str_c(pathway, collapse = ",")
  }
}

## get the drug combination, gene combination, and pathway combination
drug_data = read.csv(file='../../data/originals/drug_comb_name.csv')
drug_comb = as.character(drug_data$COMBINATION_ID)
gene_comb = rep(NA, length(drug_comb))
path_comb = rep(NA, length(drug_comb))
for (i in 1:length(drug_comb))
{
  drugs = unlist(strsplit(drug_comb[i], split = ".", fixed = TRUE))
  ind1 = which(drug == drugs[1])
  ind2 = which(drug == drugs[2])
  gene_combs = c(unlist(genelist[ind1]), unlist(genelist[ind2]))
  gene_comb[i] = str_c(gene_combs, collapse = ",")
  path_combs = paste(path_genelist[ind1], path_genelist[ind2], collapse = ",")
  path_comb[i] = gsub("[ ]", ",", path_combs)
}

# get all possible pathways
path_all = matrix()
for (i in 1:length(drug_comb))
{
  paths = unlist(strsplit(path_comb[i], split = ",", fixed = TRUE))
  paths = as.matrix(paths)
  path_all = rbind(path_all, paths)
}

path_all = unique(path_all)
path_all = sort(path_all)
# delete NAs, there are 133 possible pathways in our dataset
path_all = path_all[2:134]

## compute the score
# N is the number of genes in human
# Note that, the 'pathway.txt' was downloaded from website (I forgot the specific name)
# this file has two columns, the first one is pathway, the second one is its corresponding genes
# our purpose to use this file is to find how many genes in one specific pathway
kegg_pathway = read.table(file='../../data/originals/pathway.txt')
pathway_total = kegg_pathway$V1
hsa_total = kegg_pathway$V2
N = length(unique(hsa_total))
 
# M is the number of genes annotated to the KEGG pathway P_j
pathway_total = sapply(strsplit(as.character(pathway_total), split=':', fixed=TRUE), function(x) (x[2]))
pathway_unique = unique(pathway_total)

M = rep(NA, length(path_all))

for (i in 1:length(path_all))
{
  ind = which(pathway_unique == path_all[i])
  M[i] = length(which(pathway_total == pathway_unique[ind]))
}


# n is the number of genes in gene set G_i
n = rep(NA, length(gene_comb))
for (i in 1:length(gene_comb))
{
  n[i] = length(unlist(strsplit(gene_comb[i], split = ",", fixed = TRUE)))
}
# m is the number of genes in number of genes both in gene set G_i and in KEGG pathway P_j

m = matrix(NA, nrow = length(drug_comb), ncol = length(path_all))
for (i in 1:length(drug_comb))
{
  gene_combination = unlist(strsplit(gene_comb[i], split = ",", fixed = TRUE))
  gene_ID = paste(unlist(mget(unlist(gene_combination), org.Hs.egSYMBOL2EG, ifnotfound = NA)))
  gene_ID = paste("hsa:", gene_ID)
  gene_ID = gsub("[ ]", "", gene_ID)
  for (j in 1:length(path_all))
  {
    path_gene = as.character(hsa_total[pathway_total == path_all[j]])
    m[i, j] = length(which(gene_ID %in% path_gene == TRUE))
  }
}

# compute the score
drug_score = matrix(0, nrow = length(drug_comb), ncol = length(path_all))
for (i in 1:length(drug_comb))
  {
    for (j in 1:length(path_all))
    {
      for (k in (m[i,j]+1):(n[i]+1))
      {
        k = k-1
        temp = choose(M[j], k) * choose(N - M[j], n[i] - k) / choose(N, n[i])
        drug_score[i, j] = drug_score[i,j] + temp
      }
      drug_score[i, j] = -log10(drug_score[i,j])
    }
}


## build kernels
## correlation matrix of this 169 cell lines
drug_df = as.data.frame(t(drug_score))
cor_mat = cor(drug_df)

# dot product kernel matrix
dotprod_kernel = matrix(NA, 169, 169)

for(i in 1:169)
{
  for(j in 1:169)
  {
    l2_norm = (sqrt(t(drug_score[i,])%*%drug_score[i,])) * (sqrt(t(drug_score[j,])%*%drug_score[j,]))
    dotprod_kernel[i,j] = (drug_score[i,] %*% drug_score[j,]) / l2_norm
  }
}

# dot product angles kernel matix
angle_kernel = matrix(NA, 169, 169)
for(i in 1:169)
{
  for(j in 1:169)
  {
    if(dotprod_kernel[i, j] > 1)
    {dotprod_kernel[i, j] = 1}
    angle_kernel[i, j] = acos(dotprod_kernel[i, j]) * (180/pi)
  }
}

# rbf kernel matrix
drug_norm = matrix(NA, 169, 133)
for(i in 1:169)
{
  drug_norm[i, ]= normalize(drug_score[i, ])
}

rbf = rbfdot(sigma = 1)
rbf_kernel = matrix(NA, 169, 169)

for(i in 1:169)
{
  for(j in 1:169)
  {
    rbf_kernel[i,j] = rbf(drug_norm[i,], drug_norm[j,])
  }
}

## write files
cor_matdf = as.data.frame(round(cor_mat, 4))
write.table(cor_matdf,file='../../data/round1/kernels/corr_drug_pathway.txt',col.names=F,row.names=F,sep="\t",quote=F)

dotprod_kerneldf = as.data.frame(round(dotprod_kernel, 4))
write.table(dotprod_kerneldf,file='../../data/round1/kernels/dot_product_drug_pathway.txt',col.names=F,row.names=F,sep="\t",quote=F)

angle_kerneldf = as.data.frame(round(angle_kernel, 4))
write.table(angle_kerneldf,file='../../data/round1/kernels/angular_similarity_drug_pathway.txt',col.names=F,row.names=F,sep="\t",quote=F)

rbf_kerneldf = as.data.frame(round(rbf_kernel, 4))
write.table(rbf_kerneldf,file='../../data/round1/kernels/rbf_drug_pathway.txt',col.names=F,row.names=F,sep="\t",quote=F)









