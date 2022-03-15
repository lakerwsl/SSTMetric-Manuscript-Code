##Code for Table 3 in the manuscript
library(MASS)
library(energy)
library(class)
library(matchingR)
library(reticulate)
library(rARPACK)
library(xtable)

##Function for estimating distance Dstar
SpectralMethod <- function(X, m, n, K)
{
  d=ncol(X)
  Xbar<-matrix(0, nrow = m, ncol = d)
  SigmaSum<-matrix(0, nrow = d, ncol = d)
  for (i in 1:m)
  {
    tempX <- X[((i-1)*n+1):(i*n),]
    Xbar[i,]<-apply(tempX,2,mean)
    SigmaSum <- SigmaSum+(t(tempX) %*% tempX)/n/(n-1)
  }
  Xbarbar <- apply(Xbar,2,mean)
  Sigmabar <- (t(Xbar) %*% Xbar)/m
  Sigma<-Sigmabar*(n/(n-1)+1/(m-1))-SigmaSum/m+(Xbarbar %*% t(Xbarbar))*m/(m-1)
  Eigen<-eigs_sym(Sigma,K)
  EigenValue<-Eigen$values
  EigenValue<-EigenValue/sqrt(sum(EigenValue*EigenValue))
  Bhat<-Eigen$vectors %*% diag(sqrt(EigenValue))
  return(Bhat)
}

##Function for estimating distance Dstarstar
SpectralMethodPrac <- function(X, m, n, K)
{
  d=ncol(X)
  Xbar<-matrix(0, nrow = m, ncol = d)
  SigmaSum<-matrix(0, nrow = d, ncol = d)
  for (i in 1:m)
  {
    tempX <- X[((i-1)*n+1):(i*n),]
    Xbar[i,]<-apply(tempX,2,mean)
    SigmaSum <- SigmaSum+(t(tempX) %*% tempX)/n/(n-1)
  }
  Xbarbar <- apply(Xbar,2,mean)
  Sigmabar <- (t(Xbar) %*% Xbar)/m
  Sigma<-Sigmabar*(n/(n-1)+1/(m-1))-SigmaSum/m+(Xbarbar %*% t(Xbarbar))*m/(m-1)
  Eigen<-eigs_sym(Sigma,K)
  Bhat<-Eigen$vectors
  return(Bhat)
}

###########################################
##K=10
###########################################

## set parameters

d=100
K=10
m=1000
n=10
s=1000
lambda=4
set.seed(1)
B=matrix(0, nrow = d, ncol = K)
Bstar=matrix(0, nrow = d, ncol = K)
BstarP=matrix(0, nrow = d, ncol = K)
U=matrix(rnorm(d*d), nrow = d, ncol = d)
U=svd(U)
U=U$u
for (i in 1:K)
{
  B[,i]=U[,i]*lambda*sqrt(i/K)
  Bstar[,i]=U[,i]*sqrt(i)
  BstarP[,i]=U[,i]
}
Mstar<-Bstar %*% t(Bstar)
alpha=rep(1,K)/2
alpha[1:(K-1)]=0
Sigma=diag(rep(1,K))-alpha %*% t(alpha)

## estimate distance

set.seed(3)
mseq=c(1000,5000)
Blist<-list(4)
ii=1
for (m in mseq)
{
  Y=rbinom(m, 1, 0.5)*2-1
  Z=mvrnorm(n = m, mu=alpha, Sigma=Sigma)
  index=(Y+1)/2
  Z[index,]=-Z[index,]
  ZZ=matchingR:::reprow(Z,n)
  X=ZZ %*% t(B)
  X=X + matrix(rnorm(n*m*d), nrow = n*m, ncol = d)
  Bhat<-SpectralMethod(X, m, n, K)
  Mhat<-Bhat %*% t(Bhat)
  Blist[[ii]]<-Mhat
  Bhat<-SpectralMethodPrac(X, m, n, K)
  Mhat<-Bhat %*% t(Bhat)
  Blist[[ii+1]]<-Mhat
  ii=ii+2
}
Blist[[ii]]<-Mstar
Blist[[ii+1]]<-BstarP %*% t(BstarP)

## estimate threshold for the test

nummethod=7
times=500
Re=matrix(0,nrow = nummethod, ncol=times)
for (t in 1:times)
{
  X1=rnorm(d)
  X2=rnorm(d)
  Re[1,t]=sum((X1-X2)*(X1-X2))
  for (j in 1:(nummethod-1))
  {
    Re[j+1,t]=t((X1-X2)) %*% Blist[[j]] %*% (X1-X2)
  }
}
Threshold=rep(0,nummethod)
for(i in 1:nummethod)
{
  Threshold[i]=quantile(Re[i,],0.95)
}

## numerical experiment

ll=K/2
rr=seq(1,5,1)
powers<-matrix(0, nrow = length(rr), ncol = nummethod)
for (j in 1:length(rr))
{
  r=rr[j]
  for (t in 1:times)
  {
    X1=rnorm(d)
    X2=B %*% c(rep(r/sqrt(ll),ll),rep(0,K-ll))
    X2=X2 + rnorm(d)
    Re[1,t]=sum((X1-X2)*(X1-X2))
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=t((X1-X2)) %*% Blist[[k]] %*% (X1-X2)
    }
  }
  for(i in 1:nummethod)
  {
    powers[j,i]=mean(Re[i,]>Threshold[i])
  }
}

powers <- powers[,c(1,2,4,6,3,5,7)]

##left part of Table 3

xtable(t(powers))

###########################################
##K=50
###########################################

## set parameters

d=100
K=50
m=1000
n=10
s=1000
lambda=4
set.seed(1)
B=matrix(0, nrow = d, ncol = K)
Bstar=matrix(0, nrow = d, ncol = K)
BstarP=matrix(0, nrow = d, ncol = K)
U=matrix(rnorm(d*d), nrow = d, ncol = d)
U=svd(U)
U=U$u
for (i in 1:K)
{
  B[,i]=U[,i]*lambda*sqrt(i/K)
  Bstar[,i]=U[,i]*sqrt(i)
  BstarP[,i]=U[,i]
}
Mstar<-Bstar %*% t(Bstar)
alpha=rep(1,K)/2
alpha[1:(K-1)]=0
Sigma=diag(rep(1,K))-alpha %*% t(alpha)

## estimate distance

set.seed(3)
mseq=c(1000,5000)
Blist<-list(4)
ii=1
for (m in mseq)
{
  Y=rbinom(m, 1, 0.5)*2-1
  Z=mvrnorm(n = m, mu=alpha, Sigma=Sigma)
  index=(Y+1)/2
  Z[index,]=-Z[index,]
  ZZ=matchingR:::reprow(Z,n)
  X=ZZ %*% t(B)
  X=X + matrix(rnorm(n*m*d), nrow = n*m, ncol = d)
  Bhat<-SpectralMethod(X, m, n, K)
  Mhat<-Bhat %*% t(Bhat)
  Blist[[ii]]<-Mhat
  Bhat<-SpectralMethodPrac(X, m, n, K)
  Mhat<-Bhat %*% t(Bhat)
  Blist[[ii+1]]<-Mhat
  ii=ii+2
}
Blist[[ii]]<-Mstar
Blist[[ii+1]]<-BstarP %*% t(BstarP)

## estimate threshold for the test

nummethod=7
times=500
Re=matrix(0,nrow = nummethod, ncol=times)
for (t in 1:times)
{
  X1=rnorm(d)
  X2=rnorm(d)
  Re[1,t]=sum((X1-X2)*(X1-X2))
  for (j in 1:(nummethod-1))
  {
    Re[j+1,t]=t((X1-X2)) %*% Blist[[j]] %*% (X1-X2)
  }
}
Threshold=rep(0,nummethod)
for(i in 1:nummethod)
{
  Threshold[i]=quantile(Re[i,],0.95)
}

## numerical experiment

ll=K/2
rr=seq(1,5,1)
powers<-matrix(0, nrow = length(rr), ncol = nummethod)
for (j in 1:length(rr))
{
  r=rr[j]
  for (t in 1:times)
  {
    X1=rnorm(d)
    X2=B %*% c(rep(r/sqrt(ll),ll),rep(0,K-ll))
    X2=X2 + rnorm(d)
    Re[1,t]=sum((X1-X2)*(X1-X2))
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=t((X1-X2)) %*% Blist[[k]] %*% (X1-X2)
    }
  }
  for(i in 1:nummethod)
  {
    powers[j,i]=mean(Re[i,]>Threshold[i])
  }
}

powers <- powers[,c(1,2,4,6,3,5,7)]

##right part of Table 3

xtable(t(powers))
