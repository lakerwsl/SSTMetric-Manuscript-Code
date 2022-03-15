## Code for Table 4

library(MASS)
library(energy)
library(class)
library(matchingR)
library(reticulate)
library(rARPACK)
library(LICORS)
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

## set parameters

d=100
K=10
m=1000
n=10
s=500
lambda=2
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
alpha=c(rep(0.5/sqrt(K/2),K/2),rep(0,K/2))
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
  Blist[[ii]]<-SpectralMethod(X, m, n, K)
  Blist[[ii+1]]<-SpectralMethodPrac(X, m, n, K)
  ii=ii+2
}
Blist[[ii]]<-Bstar
Blist[[ii+1]]<-BstarP

## random start

ll=4
set.seed(5)
nummethod=7
times=500
rr=seq(0.4,1,0.2)
powerss<-matrix(0, nrow = length(rr), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(rr))
{
  r=rr[j]
  alpha=c(rep(0,K-ll),rep(r/sqrt(ll),ll))
  #alpha=c(1,rep(0,K-1))*r
  Sigma=diag(rep(1,K))-alpha %*% t(alpha)
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Truth=(Yd+1)/2+1
    
    Ree<-kmeans(Xd, centers=2, iter.max = 100, nstart = 10, algorithm="Lloyd")
    Error<-sum(Ree$cluster==Truth)/s
    Re[1,t]=min(Error,1-Error)
    
    for (k in 1:(nummethod-1))
    {
      Ree<-kmeans(Xd %*% Blist[[k]], centers=2, iter.max = 100, nstart = 10, algorithm="Lloyd")
      Error<-sum(Ree$cluster==Truth)/s
      Re[k+1,t]=min(Error,1-Error)
    }
  }
  for(i in 1:nummethod)
  {
    powerss[j,i]=mean(Re[i,])
  }
}

## left part of Table 4

powerss <- powerss[,c(1,2,4,6,3,5,7)]
xtable(t(powerss))

## perfect start

rr=seq(0.4,1,0.2)
powerss<-matrix(0, nrow = length(rr), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(rr))
{
  r=rr[j]
  alpha=c(rep(0,K-ll),rep(r/sqrt(ll),ll))
  #alpha=c(1,rep(0,K-1))*r
  Sigma=diag(rep(1,K))-alpha %*% t(alpha)
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Truth=(Yd+1)/2+1

    Center=rbind(alpha %*% t(B),-alpha %*% t(B))
    Ree<-kmeans(Xd, centers=Center, iter.max = 100, algorithm="Lloyd")
    Error<-sum(Ree$cluster==Truth)/s
    Re[1,t]=min(Error,1-Error)
    
    for (k in 1:(nummethod-1))
    {
      Center=rbind(alpha %*% t(B) %*% Blist[[k]],-alpha %*% t(B) %*% Blist[[k]])
      Ree<-kmeans(Xd %*% Blist[[k]], centers=Center, iter.max = 100, algorithm="Lloyd")
      Error<-sum(Ree$cluster==Truth)/s
      Re[k+1,t]=min(Error,1-Error)
    }
  }
  for(i in 1:nummethod)
  {
    powerss[j,i]=mean(Re[i,])
  }
}

## right part of Table 4

powerss <- powerss[,c(1,2,4,6,3,5,7)]
xtable(t(powerss))
