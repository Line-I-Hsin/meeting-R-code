FDA_sample=function(seed=123,n=20,sigma1=2,sigma2=4,t=seq(0,1,length=100),sigma3=0){
  set.seed(seed)
  a1=rnorm(n,0,sigma1)
  a2=rnorm(n,0,sigma2)
  X=matrix(0,nrow = n,ncol = length(t))
  for(i in 1:n){
    X[i,]=sqrt(2)*cos(2*pi*t)+a1[i]*sqrt(2)*sin(2*pi*t)+a2[i]*sqrt(2)*sin(3*pi*t)
  }
  e=matrix(rnorm(n*length(t),0,sigma3),nrow = n,ncol = length(t))
  Y=X+e
  ####mean
  mean_function=function(Y){
    m=col(Y)
    mean_hat=colMeans(Y)
    return(mean_hat)
  }
  ####cov
  cov_function=function(Y,mean_hat){
    n = nrow(Y)
    m = ncol(Y)
    matrix_mean_hat=matrix(mean_hat,n,m,byrow = TRUE)
    cov_hat=t(Y-matrix_mean_hat)%*%(Y-matrix_mean_hat)/(n-1)
    return(cov_hat)
  }
  ####corr
  corr_function=function(cov_hat){
    var_hat=diag(cov_hat)
    D_matrix=diag(1/sqrt(var_hat))
    corr_hat=D_matrix %*%cov_hat %*%D_matrix
    return(corr_hat)
  }
  mean_hat = mean_function(Y)
  cov_hat = cov_function(Y, mean_hat)
  corr_hat = corr_function(cov_hat)
  return(list(Y=Y,mean_hat=mean_hat,cov_hat=cov_hat,corr_hat=corr_hat,t=t))
}
result=FDA_sample()

Y=result$Y
t=result$t
mean_hat=result$mean_hat
matplot(t, t(Y), type='l', lty=1,col=rainbow(nrow(Y)),
        xlab='t', ylab='Y(t)', main='Original Data Curves')
lines(t, mean_hat, col='black',lwd=2)
library(lattice)
# 從 result 取得資料
cov_mat=result$cov_hat
corr_mat=result$corr_hat
t=result$t
#等高線圖 - 共變異數
contour(t, t, cov_mat,
        xlab = "t", ylab = "t", 
        main = "Covariance Contour Plot",
        nlevels = 20)
#等高線圖 - 相關系數
contour(t, t, corr_mat,xlab = "t", ylab = "t", main = "Correlation Contour Plot",nlevels = 20)
#3d等高線圖 - 共變異數
persp(t, t, cov_mat,theta =30, phi = 30,expand=0.5,xlab = "t", ylab = "t", zlab = "Covariance",main = "Covariance Surface Plot")
#3d等高線圖 - 相關系數
persp(t, t, corr_mat, theta=30, phi = 30,expand = 0.5,xlab = "t", ylab = "t", zlab = "Correlation",main = "Correlation Surface Plot")
####FPCA
hat_eigen=eigen(cov_mat)
var_eigen=hat_eigen$values/sum(hat_eigen$values)
caulative_persent=cumsum(var_eigen)
####
matplot(t,hat_eigen$vectors[,1:2],type="l",lty=1,xlab="t",ylab="eigenvectors")
####FPC score
n=nrow(Y)
m=ncol(Y)
matrix_mean_hat = matrix(mean_hat, nrow = n, ncol = m, byrow = TRUE)
Z_matrix=Y-matrix_mean_hat
PC_score=Z_matrix%*%hat_eigen$vectors
par(mfrow = c(1, 1))
plot(PC_score[,1],PC_score[,2],xlab="PC1",ylab="PC2",pch=16,col="blue")
points(0,0,col="red",pch=16)


result=FDA_sample(seed=123,n=20,sigma1=2,sigma2=4,t=seq(0,1,length=100),sigma3=1)
Y=result$Y
t=result$t
mean_hat=result$mean_hat
matplot(t, t(Y), type='l', lty=1,col=rainbow(nrow(Y)),
        xlab='t', ylab='Y(t)', main='Original Data Curves')
lines(t, mean_hat, col='black',lwd=2)
library(lattice)
cov_mat=result$cov_hat
corr_mat=result$corr_hat
t=result$t
#等高線圖 - 共變異數
contour(t, t, cov_mat,
        xlab = "t", ylab = "t", 
        main = "Covariance Contour Plot",
        nlevels = 20)
#等高線圖 - 相關系數
contour(t, t, corr_mat,xlab = "t", ylab = "t", main = "Correlation Contour Plot",nlevels = 20)
#3d等高線圖 - 共變異數
persp(t, t, cov_mat,theta =30, phi = 30,expand=0.5,xlab = "t", ylab = "t", zlab = "Covariance",main = "Covariance Surface Plot")
#3d等高線圖 - 相關系數
persp(t, t, corr_mat, theta=30, phi = 30,expand = 0.5,xlab = "t", ylab = "t", zlab = "Correlation",main = "Correlation Surface Plot")
####FPCA
hat_eigen=eigen(cov_mat)
var_eigen=hat_eigen$values/sum(hat_eigen$values)
caulative_persent=cumsum(var_eigen)
####
matplot(t,hat_eigen$vectors[,1:2],type="l",lty=1,xlab="t",ylab="eigenvectors")
####FPC score
n=nrow(Y)
m=ncol(Y)
matrix_mean_hat = matrix(mean_hat, nrow = n, ncol = m, byrow = TRUE)
Z_matrix=Y-matrix_mean_hat
PC_score=Z_matrix%*%hat_eigen$vectors
par(mfrow = c(1, 1))
plot(PC_score[,1],PC_score[,2],xlab="PC1",ylab="PC2",pch=16,col="blue")
points(0,0,col="red",pch=16)


result=FDA_sample(seed=123,n=20,sigma1=2,sigma2=4,t=seq(0,1,length=100),sigma3=1.5)
Y=result$Y
t=result$t
mean_hat=result$mean_hat
matplot(t, t(Y), type='l', lty=1,col=rainbow(nrow(Y)),
        xlab='t', ylab='Y(t)', main='Original Data Curves')
lines(t, mean_hat, col='black',lwd=2)
library(lattice)
cov_mat=result$cov_hat
corr_mat=result$corr_hat
#等高線圖 - 共變異數
contour(t, t, cov_mat,
        xlab = "t", ylab = "t", 
        main = "Covariance Contour Plot",
        nlevels = 20)
#等高線圖 - 相關系數
contour(t, t, corr_mat,xlab = "t", ylab = "t", main = "Correlation Contour Plot",nlevels = 20)
#3d等高線圖 - 共變異數
persp(t, t, cov_mat,theta =30, phi = 30,expand=0.5,xlab = "t", ylab = "t", zlab = "Covariance",main = "Covariance Surface Plot")
#3d等高線圖 - 相關系數
persp(t, t, corr_mat, theta=30, phi = 30,expand = 0.5,xlab = "t", ylab = "t", zlab = "Correlation",main = "Correlation Surface Plot")
####FPCA
hat_eigen=eigen(cov_mat)
var_eigen=hat_eigen$values/sum(hat_eigen$values)
caulative_persent=cumsum(var_eigen)
####
matplot(t,hat_eigen$vectors[,1:4],type="l",lty=1,xlab="t",ylab="eigenvectors")
####FPC score
n=nrow(Y)
m=ncol(Y)
matrix_mean_hat = matrix(mean_hat, nrow = n, ncol = m, byrow = TRUE)
Z_matrix=Y-matrix_mean_hat
PC_score=Z_matrix%*%hat_eigen$vectors
par(mfrow = c(2, 3))
plot(PC_score[,1],PC_score[,2],xlab="PC1",ylab="PC2",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,3],xlab="PC1",ylab="PC3",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,4],xlab="PC1",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,3],xlab="PC2",ylab="PC3",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,4],xlab="PC2",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,3],PC_score[,4],xlab="PC3",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)



result=FDA_sample(seed=123,n=20,sigma1=2,sigma2=4,t=seq(0,1,length=100),sigma3=2)
Y=result$Y
t=result$t
mean_hat=result$mean_hat
matplot(t, t(Y), type='l', lty=1,col=rainbow(nrow(Y)),
        xlab='t', ylab='Y(t)', main='Original Data Curves')
lines(t, mean_hat, col='black',lwd=2)
library(lattice)
cov_mat=result$cov_hat
corr_mat=result$corr_hat
#等高線圖 - 共變異數
contour(t, t, cov_mat,
        xlab = "t", ylab = "t", 
        main = "Covariance Contour Plot",
        nlevels = 20)
#等高線圖 - 相關系數
contour(t, t, corr_mat,xlab = "t", ylab = "t", main = "Correlation Contour Plot",nlevels = 20)
#3d等高線圖 - 共變異數
persp(t, t, cov_mat,theta =30, phi = 30,expand=0.5,xlab = "t", ylab = "t", zlab = "Covariance",main = "Covariance Surface Plot")
#3d等高線圖 - 相關系數
persp(t, t, corr_mat, theta=30, phi = 30,expand = 0.5,xlab = "t", ylab = "t", zlab = "Correlation",main = "Correlation Surface Plot")
####FPCA
hat_eigen=eigen(cov_mat)
var_eigen=hat_eigen$values/sum(hat_eigen$values)
caulative_persent=cumsum(var_eigen)
####
matplot(t,hat_eigen$vectors[,1:8],type="l",lty=1,xlab="t",ylab="eigenvectors")
####FPC score
n=nrow(Y)
m=ncol(Y)
matrix_mean_hat = matrix(mean_hat, nrow = n, ncol = m, byrow = TRUE)
Z_matrix=Y-matrix_mean_hat
PC_score=Z_matrix%*%hat_eigen$vectors
par(mfrow = c(2, 4))
plot(PC_score[,1],PC_score[,2],xlab="PC1",ylab="PC2",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,3],xlab="PC1",ylab="PC3",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,4],xlab="PC1",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,5],xlab="PC1",ylab="PC5",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,6],xlab="PC1",ylab="PC6",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,7],xlab="PC1",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,1],PC_score[,8],xlab="PC1",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,3],xlab="PC2",ylab="PC3",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,4],xlab="PC2",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,5],xlab="PC2",ylab="PC5",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,6],xlab="PC2",ylab="PC6",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,7],xlab="PC2",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,2],PC_score[,8],xlab="PC2",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
par(mfrow = c(2, 3))
plot(PC_score[,3],PC_score[,4],xlab="PC3",ylab="PC4",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,3],PC_score[,5],xlab="PC3",ylab="PC5",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,3],PC_score[,6],xlab="PC3",ylab="PC6",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,3],PC_score[,7],xlab="PC3",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,3],PC_score[,8],xlab="PC3",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
par(mfrow = c(2, 2))
plot(PC_score[,4],PC_score[,5],xlab="PC4",ylab="PC5",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,4],PC_score[,6],xlab="PC4",ylab="PC6",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,4],PC_score[,7],xlab="PC4",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,4],PC_score[,8],xlab="PC4",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
par(mfrow = c(2, 3))
plot(PC_score[,5],PC_score[,6],xlab="PC5",ylab="PC6",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,5],PC_score[,7],xlab="PC5",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,5],PC_score[,8],xlab="PC5",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,6],PC_score[,7],xlab="PC6",ylab="PC7",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,6],PC_score[,8],xlab="PC6",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
plot(PC_score[,7],PC_score[,8],xlab="PC7",ylab="PC8",pch=16,col="blue")
points(0,0,col="red",pch=16)
