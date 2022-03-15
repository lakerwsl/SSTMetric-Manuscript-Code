## Code for Figure 4

library(MASS)
library(energy)
library(class)
library(matchingR)
library(reticulate)
library(rARPACK)

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
s=1000
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
alpha=c(rep(0.9/sqrt(K/2),K/2),rep(0,K/2))
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

## numerical experiments for different sample size

nummethod=7
times=100
ss=seq(500,5000,500)
sss=500
powers<-matrix(0, nrow = length(ss), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(ss))
{
  s=ss[j]
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Yd=(Yd+1)/2
    
    Yt=rbinom(sss, 1, 0.5)*2-1
    Zt=mvrnorm(n = sss, mu=alpha, Sigma=Sigma)
    index=as.logical((Yt+1)/2)
    Zt[index,]=-Zt[index,]
    Xt=Zt %*% t(B)
    Xt=Xt + matrix(rnorm(sss*d), nrow = sss, ncol = d)
    Yt=(Yt+1)/2
    Truth=Yt
    
    beta=6
    Neighbour1=ceiling(s^(beta/(beta+d)))*2
    Neighbour2=ceiling(s^(beta/(beta+K)))*2
    
    Ree<-class:::knn(train = Xd, test = Xt, cl = factor(Yd), k = Neighbour1)
    Ree<-as.numeric(as.character(Ree))
    Re[1,t]<-sum(Ree==Truth)/sss
    
    for (k in 1:(nummethod-1))
    {
      Ree<-class:::knn(train = Xd %*% Blist[[k]], test = Xt %*% Blist[[k]], cl = factor(Yd), k = Neighbour2)
      Ree<-as.numeric(as.character(Ree))
      Re[k+1,t]<-sum(Ree==Truth)/sss
    }
  }
  for(i in 1:nummethod)
  {
    powers[j,i]=mean(Re[i,])
  }
}

powerss=1-powers
PlotData<-data.frame(Error=c(powerss[,1],powerss[,2],powerss[,4],powerss[,6],powerss[,3],powerss[,5],powerss[,7]),s=rep(ss,7),Method=c(rep("Euclidean distance",10),rep("Dstar (m=1000)",10),rep("Dstar (m=5000)",10),rep("Dstar (true)",10),rep("Dstarstar (m=1000)",10),rep("Dstarstar (m=5000)",10),rep("Dstarstar (true)",10)))

ggplot(data = PlotData, aes(x=s, y=Error, group=Method, color=Method)) +
  geom_line() + theme(legend.position="bottom")

pp1<-ggplot(data = PlotData, aes(x=s, y=Error, group=Method, color=Method)) +
  geom_line() + theme(legend.position="bottom") + ylim(0.15, 0.5)

## numerical experiments for different r

rr=seq(0.1,1,0.1)
s=2000
sss=500
powers2<-matrix(0, nrow = length(rr), ncol = nummethod)
Re=matrix(0,nrow = nummethod, ncol=times)
for (j in 1:length(rr))
{
  r=rr[j]
  alpha=c(rep(r/sqrt(K/2),K/2),rep(0,K/2))
  Sigma=diag(rep(1,K))-alpha %*% t(alpha)
  for (t in 1:times)
  {
    Yd=rbinom(s, 1, 0.5)*2-1
    Zd=mvrnorm(n = s, mu=alpha, Sigma=Sigma)
    index=as.logical((Yd+1)/2)
    Zd[index,]=-Zd[index,]
    Xd=Zd %*% t(B)
    Xd=Xd + matrix(rnorm(s*d), nrow = s, ncol = d)
    Yd=(Yd+1)/2
    
    Yt=rbinom(sss, 1, 0.5)*2-1
    Zt=mvrnorm(n = sss, mu=alpha, Sigma=Sigma)
    index=as.logical((Yt+1)/2)
    Zt[index,]=-Zt[index,]
    Xt=Zt %*% t(B)
    Xt=Xt + matrix(rnorm(sss*d), nrow = sss, ncol = d)
    Yt=(Yt+1)/2
    Truth=Yt
    
    beta=6
    Neighbour1=ceiling(s^(beta/(beta+d)))*2
    Neighbour2=ceiling(s^(beta/(beta+K)))*2
    
    Ree<-class:::knn(train = Xd, test = Xt, cl = factor(Yd), k = Neighbour1)
    Ree<-as.numeric(as.character(Ree))
    Re[1,t]<-sum(Ree==Truth)/sss
    
    for (k in 1:(nummethod-1))
    {
      Ree<-class:::knn(train = Xd %*% Blist[[k]], test = Xt %*% Blist[[k]], cl = factor(Yd), k = Neighbour2)
      Ree<-as.numeric(as.character(Ree))
      Re[k+1,t]<-sum(Ree==Truth)/sss
    }
  }
  for(i in 1:nummethod)
  {
    powers2[j,i]=mean(Re[i,])
  }
}

powerss=1-powers2
PlotData<-data.frame(Error=c(powerss[,1],powerss[,2],powerss[,4],powerss[,6],powerss[,3],powerss[,5],powerss[,7]),r=rep(rr,7),Method=c(rep("Euclidean distance",10),rep("Dstar (m=1000)",10),rep("Dstar (m=5000)",10),rep("Dstar (true)",10),rep("Dstarstar (m=1000)",10),rep("Dstarstar (m=5000)",10),rep("Dstarstar (true)",10)))

ggplot(data = PlotData, aes(x=r, y=Error, group=Method, color=Method)) +
  geom_line() + theme(legend.position="bottom")

legendtex <- list()
legendtex[[1]] <- TeX(r'($D^{*} (m=1000)$)')
legendtex[[2]] <- TeX(r'($D^{*} (m=5000)$)')
legendtex[[3]] <- TeX(r'($D^{*} (true)$)')
legendtex[[4]] <- TeX(r'($D^{**} (m=1000)$)')
legendtex[[5]] <- TeX(r'($D^{**} (m=5000)$)')
legendtex[[6]] <- TeX(r'($D^{**} (true)$)')
legendtex[[7]] <- TeX(r'(Euclidean)')

ggplot(data = PlotData, aes(x=r, y=Error, group=Method, color=Method)) +
  geom_line() + theme(legend.position="bottom") +
  scale_color_discrete(labels=legendtex) 

pp2<-ggplot(data = PlotData, aes(x=r, y=Error, group=Method, color=Method)) +
  geom_line() + theme(legend.position="bottom") + ylim(0.15, 0.5) +
  scale_color_discrete(labels=legendtex) 

## combine figures -- Figure 4

legend <- get_legend(pp2)
pp1 <- pp1 + theme(legend.position="none")
pp2 <- pp2 + theme(legend.position="none")
grid.arrange(pp1, pp2, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.3))
