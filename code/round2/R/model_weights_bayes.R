## Homework 8 

# 9.2.a
data=read.table('azdiabetes.dat',header=T)

# Input
data=data[,-ncol(data)]
y=as.matrix(data[,2])
X=as.matrix(data[,-2])

# Parameters
g=length(y)
nu0=2
s20=1

# N of independent samples
S=1000

# Code from page 159
n=dim(X)[1]
p=dim(X)[2]
Hg= (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
SSRg= t(y)%*%( diag(1,nrow=n) - Hg )%*%y

s2=1/rgamma(S,(nu0+n)/2,(nu0*s20+SSRg)/2)

Vb=g*solve(t(X)%*%X)/(g+1)
Eb=Vb%*%t(X)%*%y

E=matrix(rnorm(S*p,0,sqrt(s2)),S,p)
beta=t(t(E%*%chol(Vb))+c(Eb))

# 9.2.b
betas_matrix=c()
# Function to compute marginal probability
lpy.X=function(y,X,g=length(y),nu0=2,s20=try(summary(lm(y~-1+X))$sigma^2,silent=TRUE))
{
    n=dim(X)[1]
    p=dim(X)[2]
    if(p==0){Hg=0;s20=mean(y^2)}
    if(p>0) {Hg= (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)}
    SSRg= t(y)%*%( diag(1,nrow=n) - Hg )%*%y 

    S=1
    s2=1/rgamma(S,(nu0+n)/2,(nu0*s20+SSRg)/2)

    Vb=g*solve(t(X)%*%X)/(g+1)
    Eb=Vb%*%t(X)%*%y

    E=matrix(rnorm(S*p,0,sqrt(s2)),S,p)
    betas=t(t(E%*%chol(Vb))+c(Eb))
    betas_matrix=rbind(betas_matrix,betas)

    -0.5*(n*log(pi)+p*log(1+g)+(nu0+n)*log(nu0*s20+SSRg)-nu0*log(nu0*s20))+lgamma((nu0+n)/2)-lgamma(nu0/2)

}

# MCMC starting values
z=rep(1,dim(X)[2])
lpy.c=lpy.X(y,X[,z==1,drop=FALSE])
S=1000
Z=matrix(NA,S,dim(X)[2])

# Gibbs sampler
for(s in 1:S)
{
    for(j in sample(1:dim(X)[2]))
    {
        zp=z;zp[j]=1-zp[j]
        lpy.p=lpy.X(y,X[,zp==1,drop=FALSE])
        r=(lpy.p-lpy.c)*(-1)^(zp[j]==0)
        z[j]=rbinom(1,1,1/(1+exp(-r)))
        if(z[j]==zp[j]) {lpy.c=lpy.p}
    }
    Z[s,]=z
}

## 9.3.a
data=read.table('crime.dat',header=T)

# Input
y=as.matrix(data[,1])
X=as.matrix(data[,-1])

# Parameters
g=length(y)
nu0=2
s20=1

# N of independent samples
S=1000

# Code from page 159
n=dim(X)[1]
p=dim(X)[2]
Hg= (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
SSRg= t(y)%*%( diag(1,nrow=n) - Hg )%*%y

s2=1/rgamma(S,(nu0+n)/2,(nu0*s20+SSRg)/2)

Vb=g*solve(t(X)%*%X)/(g+1)
Eb=Vb%*%t(X)%*%y

E=matrix(rnorm(S*p,0,sqrt(s2)),S,p)
beta=t(t(E%*%chol(Vb))+c(Eb))

# 9.3.b
train=data[1:(nrow(data)/2),]
test=data[(1+nrow(data)/2):nrow(data),]
y_train=as.matrix(train[,1])
X_train=as.matrix(train[,-1])
y_test=as.matrix(test[,1])
X_test=as.matrix(test[,-1])

fit=lm(y~M+So+Ed+Po1+Po2+LF+M.F+Pop+NW+U1+U2+GDP+Ineq+Prob+Time, data=train)
prediction=predict(fit,data.frame(X_test))
plot(test[,1],prediction)

ols_mse=0
for(i in 1:nrow(test))
{
    ols_mse=ols_mse+(test[i,1]-prediction[i])^2
}
ols_mse=ols_mse/nrow(test)


# Bayesian prediction
# Parameters
g=length(y_train)
nu0=2
s20=1

# N of independent samples
S=1000

# Code from page 159
n=dim(X_train)[1]
p=dim(X_train)[2]
Hg= (g/(g+1)) * X_train%*%solve(t(X_train)%*%X_train)%*%t(X_train)
SSRg= t(y_train)%*%( diag(1,nrow=n) - Hg )%*%y_train

s2=1/rgamma(S,(nu0+n)/2,(nu0*s20+SSRg)/2)

Vb=g*solve(t(X_train)%*%X_train)/(g+1)
Eb=Vb%*%t(X_train)%*%y_train

E=matrix(rnorm(S*p,0,sqrt(s2)),S,p)
beta=t(t(E%*%chol(Vb))+c(Eb))

beta_vector=vector(length=ncol(beta))
for(i in 1:ncol(beta))
{
    beta_vector[i]=mean(beta[,i])
}

bayes_prediction=X_test%*%beta_vector

plot(test[,1],bayes_prediction)

bayes_mse=0
for(i in 1:nrow(test))
{
    bayes_mse=bayes_mse+(test[i,1]-bayes_prediction[i])^2
}
bayes_mse=bayes_mse/nrow(test)

# 9.3.c
Niter=100
ols_mse_vector=vector(length=Niter)
bayes_mse_vector=vector(length=Niter)

smp_size=floor(0.5*nrow(data))
set.seed(123)

for( j in 1:Niter)
{
    # Generate random sets
    train_ind=sample(seq_len(nrow(data)), size = smp_size)
    train=data[train_ind, ]
    test=data[-train_ind, ]
    y_train=as.matrix(train[,1])
    X_train=as.matrix(train[,-1])
    y_test=as.matrix(test[,1])
    X_test=as.matrix(test[,-1])
    ## OLS
    fit=lm(y~M+So+Ed+Po1+Po2+LF+M.F+Pop+NW+U1+U2+GDP+Ineq+Prob+Time, data=train)
    prediction=predict(fit,data.frame(X_test))
    ols_mse=0
    for(i in 1:nrow(test))
    {
        ols_mse=ols_mse+(test[i,1]-prediction[i])^2
    }
    ols_mse=ols_mse/nrow(test)
    ## Bayes
    g=length(y_train)
    nu0=2
    s20=1
    # N of independent samples
    S=1000
    # Code from page 159
    n=dim(X_train)[1]
    p=dim(X_train)[2]
    Hg= (g/(g+1)) * X_train%*%solve(t(X_train)%*%X_train)%*%t(X_train)
    SSRg= t(y_train)%*%( diag(1,nrow=n) - Hg )%*%y_train
    s2=1/rgamma(S,(nu0+n)/2,(nu0*s20+SSRg)/2)
    Vb=g*solve(t(X_train)%*%X_train)/(g+1)
    Eb=Vb%*%t(X_train)%*%y_train
    E=matrix(rnorm(S*p,0,sqrt(s2)),S,p)
    beta=t(t(E%*%chol(Vb))+c(Eb))
    beta_vector=vector(length=ncol(beta))
    for(i in 1:ncol(beta))
    {
        beta_vector[i]=mean(beta[,i])
    }
    bayes_prediction=X_test%*%beta_vector
    bayes_mse=0
    for(i in 1:nrow(test))
    {
        bayes_mse=bayes_mse+(test[i,1]-bayes_prediction[i])^2
    }
    bayes_mse=bayes_mse/nrow(test)
    
    ols_mse_vector[j]=ols_mse
    bayes_mse_vector[j]=bayes_mse
}

diff_ols_bayes=ols_mse_vector-bayes_mse_vector
plot(density(ols_mse_vector-bayes_mse_vector))
