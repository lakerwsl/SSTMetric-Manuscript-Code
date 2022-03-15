## Code for right part of Table 5

library(keras)
library(EBImage)
library(abind)
library(foreach)
library(doParallel)
library(doRNG)
library(xtable)

##Find data
fashion_mnist <- dataset_fashion_mnist()
xx_train <- fashion_mnist$train$x
yy_train <- fashion_mnist$train$y
xx_test <- fashion_mnist$test$x
yy_test <- fashion_mnist$test$y

rotate <- function(x) t(apply(x, 2, rev))
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
  Eigen<-eigen(Sigma)
  EigenValue<-Eigen$values[1:K]
  EigenValue<-EigenValue/sqrt(sum(EigenValue*EigenValue))
  Bhat<-Eigen$vectors[,1:K] %*% diag(sqrt(EigenValue))
  #Bhat<-Eigen$vectors[,1:K]
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
  Eigen<-eigen(Sigma)
  EigenValue<-Eigen$values[1:K]
  EigenValue<-EigenValue/sqrt(sum(EigenValue*EigenValue))
  #Bhat<-Eigen$vectors[,1:K] %*% diag(sqrt(EigenValue))
  Bhat<-Eigen$vectors[,1:K]
  return(Bhat)
}

##Function for augmentation
DataAugmentation <- function(X, drift)
{
  SS=dim(X)
  numimg=SS[1]
  DataX=matrix(0,nrow = numimg*5, ncol=SS[2]*SS[3])
  for (i in 1:numimg)
  {
    Timg=X[i,,]
    start=(i-1)*5
    DataX[start+1,]=c(Timg)/255
    DataX[start+2,]=c(rbind(Timg[(SS[2]-drift+1):SS[2],],Timg[1:(SS[2]-drift),]))/255
    DataX[start+3,]=c(rbind(Timg[(drift+1):SS[2],],Timg[1:drift,]))/255
    DataX[start+4,]=c(cbind(Timg[,(SS[3]-drift+1):SS[3]],Timg[,1:(SS[3]-drift)]))/255
    DataX[start+5,]=c(cbind(Timg[,(drift+1):SS[3]],Timg[,1:drift]))/255
  }
  return(DataX)
}

DataTransformation <- function(X)
{
  SS=dim(X)
  numimg=SS[1]
  DataX=matrix(0,nrow = numimg, ncol=SS[2]*SS[3])
  for (i in 1:numimg)
  {
    Timg=X[i,,]
    DataX[i,]=c(Timg)/255
  }
  return(DataX)
}

SelfTrainSize=10000
TestSize=1000

times=10
cl <- makeCluster(times) #not to overload your computer
registerDoParallel(cl)

set.seed(7)
dim=80
TS <- c(1000,2000,5000)
Ree <- matrix(0, nrow = times, ncol = 3)
Errorr <- matrix(0, nrow = 3, ncol = length(TS))
for (j in 1:length(TS))
{
  TrainSize=TS[j]
  Ree<-foreach(t=1:times, .combine=rbind, .packages='class') %dorng% {
    Reee=rep(0,3)
    SelfTrainIndex=sample(1:60000,SelfTrainSize+TrainSize)
    TrainIndex=SelfTrainIndex[1:TrainSize]
    SelfTrainIndex=SelfTrainIndex[-(1:TrainSize)]
    TestIndex=sample(1:10000,TestSize)
    
    SelfTrainX=DataAugmentation(xx_train[SelfTrainIndex,,],2)
    TrainX=DataTransformation(xx_train[TrainIndex,,])
    TestX=DataTransformation(xx_test[TestIndex,,])
    TrainLabel=yy_train[TrainIndex]
    TestLabel=yy_test[TestIndex]
    
    Bhar1<-SpectralMethod(SelfTrainX, SelfTrainSize, 5, dim)
    Bhar2<-SpectralMethodPrac(SelfTrainX, SelfTrainSize, 5, dim)
    
    ##Euclidean distance
    Re<-knn(train = TrainX, test = TestX, cl = factor(TrainLabel), k = 3)
    Truth=factor(TestLabel)
    Reee[1]<-sum(Re==Truth)/length(TestLabel)
    
    ##D star
    Re<-knn(train = TrainX %*% Bhar1, test = TestX %*% Bhar1, cl = factor(TrainLabel), k = 5)
    Truth=factor(TestLabel)
    Reee[2]<-sum(Re==Truth)/length(TestLabel)
    
    ##D star star
    Re<-knn(train = TrainX %*% Bhar2, test = TestX %*% Bhar2, cl = factor(TrainLabel), k = 5)
    Truth=factor(TestLabel)
    Reee[3]<-sum(Re==Truth)/length(TestLabel)
    
    Reee
  }
  Errorr[1,j]<-mean(Ree[,1])
  Errorr[2,j]<-mean(Ree[,2])
  Errorr[3,j]<-mean(Ree[,3])
}

stopCluster(cl)

##Results: each row corresponds to sample size and each column corresponds to each method
##left part of Table 5

xtable(t(1-Errorr),digits =3)





