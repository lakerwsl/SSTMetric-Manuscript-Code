## Code for Figure 3
library(MASS)
library(energy)
library(class)
library(matchingR)
library(reticulate)
library(rARPACK)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(grid)
library(latex2exp)

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

##Function for energy distance test
Twosample <- function(X,Y)
{
  indexplus=as.logical((Y+1)/2)
  indexminus=as.logical((-Y+1)/2)
  nplus=sum(indexplus)
  nminus=sum(indexminus)
  Xplus=X[indexplus,]
  Xminus=X[indexminus,]
  
  Xplusbar<-apply(Xplus,2,mean)
  Splus<-sum(Xplusbar*Xplusbar)*nplus/(nplus-1)
  TempS=0
  for (i in 1:nplus)
  {
    TempS <- TempS+sum(Xplus[i,]*Xplus[i,])/nplus
  }
  Splus<-Splus-TempS/(nplus-1)
  
  Xminusbar<-apply(Xminus,2,mean)
  Sminus<-sum(Xminusbar*Xminusbar)*nminus/(nminus-1)
  TempS=0
  for (i in 1:nminus)
  {
    TempS <- TempS+sum(Xminus[i,]*Xminus[i,])/nminus
  }
  Sminus<-Sminus-TempS/(nminus-1)
  
  Scross<-sum(Xplusbar*Xminusbar)*2
  
  Stat<-Splus+Sminus-Scross
  return(Stat)
}

## set parameters

d=100
K=10
m=1000
n=10
s=500
lambda=1
set.seed(1)
B=matrix(0, nrow = d, ncol = K)
Bstar=matrix(0, nrow = d, ncol = K)
U=matrix(rnorm(d*d), nrow = d, ncol = d)
BstarP=matrix(0, nrow = d, ncol = K)
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


nummethod=7
times=500
Re=matrix(0,nrow = nummethod, ncol=times)
ll=6

## numerical experiments for different r

rr=seq(0,0.5,0.05)
powers<-matrix(0, nrow = length(rr), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(rr))
{
  r=rr[j]
  alpha=c(rep(r/sqrt(ll),ll),rep(0,K-ll))
  Sigma=diag(rep(1,K))-alpha %*% t(alpha)
  
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=rep(0,K), Sigma=Sigma)
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Re[1,t]=Twosample(Xd,Yd)
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=Twosample(Xd %*% Blist[[k]],Yd)
    }
  }
  Threshold=rep(0,nummethod)
  for(i in 1:nummethod)
  {
    Threshold[i]=quantile(Re[i,],0.95)
  }
  
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Re[1,t]=Twosample(Xd,Yd)
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=Twosample(Xd %*% Blist[[k]],Yd)
    }
  }
  for(i in 1:nummethod)
  {
    powers[j,i]=mean(Re[i,]>Threshold[i])
  }
}

powers

PlotData<-data.frame(Power=c(powers[,1],powers[,2],powers[,4],powers[,6],powers[,3],powers[,5],powers[,7]),r=rep(rr,7),Method=c(rep("Euclidean distance",11),rep("Dstar (m=1000)",11),rep("Dstar (m=5000)",11),rep("Dstar (true)",11),rep("Dstarstar (m=1000)",11),rep("Dstarstar (m=5000)",11),rep("Dstarstar (true)",11)))
ggplot(data = PlotData, aes(x=r, y=Power, group=Method, color=Method)) +
  geom_line() + ylim(0, 1)

legendtex <- list()
legendtex[[1]] <- TeX(r'($D^{*} (m=1000)$)')
legendtex[[2]] <- TeX(r'($D^{*} (m=5000)$)')
legendtex[[3]] <- TeX(r'($D^{*} (true)$)')
legendtex[[4]] <- TeX(r'($D^{**} (m=1000)$)')
legendtex[[5]] <- TeX(r'($D^{**} (m=5000)$)')
legendtex[[6]] <- TeX(r'($D^{**} (true)$)')
legendtex[[7]] <- TeX(r'(Euclidean)')

ggplot(data = PlotData, aes(x=r, y=Power, group=Method, color=Method)) +
  geom_line() + ylim(0, 1) +
  scale_color_discrete(labels=legendtex) 

pp1<-ggplot(data = PlotData, aes(x=r, y=Power, group=Method, color=Method)) +
  geom_line() +
  ylim(0, 1) + 
  theme(legend.position="bottom") +
  scale_color_discrete(labels=legendtex) 

## numerical experiments for different lambda

lambdas=seq(0.5,5,0.5)
powerss<-matrix(0, nrow = length(lambdas), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(lambdas))
{
  lambda=lambdas[j]
  for (i in 1:K)
  {
    B[,i]=U[,i]*lambda*sqrt(i/K)
  }
  alpha=c(rep(1/sqrt(ll),ll),rep(0,K-ll))*0.3/lambda
  Sigma=diag(rep(1,K))-alpha %*% t(alpha)
  
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=rep(0,K), Sigma=Sigma)
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Re[1,t]=Twosample(Xd,Yd)
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=Twosample(Xd %*% Blist[[k]],Yd)
    }
  }
  Threshold=rep(0,nummethod)
  for(i in 1:nummethod)
  {
    Threshold[i]=quantile(Re[i,],0.95)
  }
  
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Re[1,t]=Twosample(Xd,Yd)
    for (k in 1:(nummethod-1))
    {
      Re[k+1,t]=Twosample(Xd %*% Blist[[k]],Yd)
    }
  }
  for(i in 1:nummethod)
  {
    powerss[j,i]=mean(Re[i,]>Threshold[i])
  }
}

powerss
PlotData<-data.frame(Power=c(powerss[,1],powerss[,2],powerss[,4],powerss[,6],powerss[,3],powerss[,5],powerss[,7]),lambda=rep(lambdas,7),Method=c(rep("Euclidean distance",10),rep("Dstar (m=1000)",10),rep("Dstar (m=5000)",10),rep("Dstar (true)",10),rep("Dstarstar (m=1000)",10),rep("Dstarstar (m=5000)",10),rep("Dstarstar (true)",10)))
ggplot(data = PlotData, aes(x=lambda, y=Power, group=Method, color=Method)) +
  geom_line() +
  xlab(TeX(r'($\lambda$)')) +
  ylim(0, 1) 
  

pp2<-ggplot(data = PlotData, aes(x=lambda, y=Power, group=Method, color=Method)) +
  geom_line() +
  ylim(0, 1) +
  theme(legend.position="bottom") +
  xlab(TeX(r'($\lambda$)'))

## combine figures -- Figure 3

legend <- get_legend(pp1)
pp1 <- pp1 + theme(legend.position="none")
pp2 <- pp2 + theme(legend.position="none")
grid.arrange(pp1, pp2, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.3))
