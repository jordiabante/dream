library(SpatialEpi)
library(kernlab)
library(BioPhysConnectoR)

cnv_gene = read.csv('cnv_gene.csv')
cnv_segment = read.csv('cnv_segment.csv')
therapy = read.csv('therapy.csv')

attach(cnv_segment)
cell_line_org = unique(cell_line_name)
cell_line = cell_line_org[order(cell_line_org)]

combination = cbind(minorCN, totalCN)
combination = unique(combination)
combination = mat.sort(combination, 1, decreasing = FALSE)
comb = matrix(NA, nrow(combination), 2)
comb[,1] = c(rep(0, 15), rep(1, 13), rep(2, 11), rep(3, 9), rep(4, 7), rep(5, 2), rep(6,1))
comb[,2] = c(seq(0,14), seq(2,14), seq(4, 14), seq(6, 14), seq(8,14), c(10,14), 12)

cn_features = matrix(NA, 85, nrow(comb))

for (i in 1:85)
{
  CN = cbind(minorCN[cell_line_name == cell_line[i]], totalCN[cell_line_name == cell_line[i]])
  for (j in 1:nrow(comb))
  {
    cn_features[i, j] = length(intersect(which(CN[,1] == comb[j, 1]), which(CN[,2] == comb[j, 2])))
  }
  #cn_features[i,] = cn_features[i,]/sum(cn_features[i,])
  
}


## correlation matrix of this 85 cell lines
cn_df = as.data.frame(t(cn_features))
cor_mat = cor(cn_df)

# dot product kernel matrix
dotprod_kernel = matrix(NA, 85, 85)

for(i in 1:85)
{
  for(j in 1:85)
  {
    l2_norm = (sqrt(t(cn_features[i,])%*%cn_features[i,])) * (sqrt(t(cn_features[j,])%*%cn_features[j,]))
    dotprod_kernel[i,j] = (cn_features[i,] %*% cn_features[j,]) / l2_norm
  }
}

# rbf kernel matrix
cn_norm = matrix(NA, 85, 58)
for(i in 1:85)
{
  cn_norm[i, ]= normalize(cn_features[i, ])
}

rbf = rbfdot(sigma = 1)
rbf_kernel = matrix(NA, 85, 85)

for(i in 1:85)
{
  for(j in 1:85)
  {
    rbf_kernel[i,j] = rbf(cn_norm[i,], cn_norm[j,])
  }
}

# write files
cor_matdf = as.data.frame(cor_mat)
colnames(cor_matdf) = cell_line
rownames(cor_matdf) = cell_line
write.csv(cor_matdf, 'cor_mat.csv')

dotprod_kerneldf = as.data.frame(dotprod_kernel)
colnames(dotprod_kerneldf) = cell_line
rownames(dotprod_kerneldf) = cell_line
write.csv(dotprod_kerneldf, 'dotprod_kernel.csv')

rbf_kerneldf = as.data.frame(rbf_kernel)
colnames(rbf_kerneldf) = cell_line
rownames(rbf_kerneldf) = cell_line
write.csv(rbf_kerneldf, 'rbf_kernel.csv')




