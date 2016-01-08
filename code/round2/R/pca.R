# Read in data
X=read.table('../../data/round2/predictor_matrix.txt',header=T)
X=as.matrix(X)
X_c=scale(X,center=T,scale=F)

# PCA
pc=princomp(X_c)
pc_1= X_c %*% pc$loadings[,1]

# SVD
n_pred=10
ss=svd(X_c)
dd=round(ss$d^2/sum(ss$d^2),2)
U=ss$u[,1:n_pred]
D=diag(ss$d[1:n_pred])
V=ss$v[1:n_pred,]

X_filt=U%*%D%*%t(V)
