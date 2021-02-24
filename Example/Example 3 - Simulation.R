# Simulation with zero nuisance parameters 
library(latex2exp)
## NULL Setting
N<-1000
n<-400 # 100, 200, 300, 400, 500

ME3<-matrix(0, N, 3)
A3<-matrix(0, N, 3)
B3<-matrix(0, N, 3)
ALPHA3<-matrix(0, N, 3)

# Simulation 
for (f in 1:N)
{
  ## True parameters
  beta1<-c(1.2, 0, 0.5)
  beta2<-c(1.2, 0, -0.8)
  a<-c(0.5, 0, 0.7)
  b<-c(0.6, 0, -0.5)
  c<-c(-1, 0, 1)
  sigma1<-c(0.6, 0.4, 0.5)
  sigma2<-c(0.5, 0.7, 0.6)
  alpha<-c(0.2, 0.5, 0.3)
  ## Sample 
  Y<-NULL
  M<-NULL
  n1<-floor(n*alpha[1])
  n2<-floor(n*alpha[2])
  n3<-n-n1-n2
  X<-rnorm(n)
  M[1:n1]<-beta1[1]+a[1]*X[1:n1]+rnorm(n1)*sigma1[1]
  Y[1:n1]<-beta2[1]+c[1]*X[1:n1]+b[1]*M[1:n1]+rnorm(n1)*sigma2[1]
  M[(1+n1):(n1+n2)]<-beta1[2]+a[2]*X[(1+n1):(n1+n2)]+rnorm(n2)*sigma1[2]
  Y[(1+n1):(n1+n2)]<-beta2[2]+c[2]*X[(1+n1):(n1+n2)]+b[2]*M[(1+n1):(n1+n2)]+rnorm(n2)*sigma2[2]
  M[(1+n1+n2):n]<-beta1[3]+a[3]*X[(1+n1+n2):n]+rnorm(n3)*sigma1[3]
  Y[(1+n1+n2):n]<-beta2[3]+c[3]*X[(1+n1+n2):n]+b[3]*M[(1+n1+n2):n]+rnorm(n3)*sigma2[3]
  U<-rep(1, n)
  V<-matrix(c(U, X), ncol=2)
  Z<-matrix(c(U, X, M), ncol=3)
  
  # EM Estimation for K=3
  ## Initialization parameters
  K<-3
  beta1<-runif(K)
  beta2<-runif(K)
  a<-runif(K)
  b<-runif(K)
  c<-runif(K)
  d1<-runif(1, 0.3, 0.4)
  d2<-runif(1, 0.3, 0.4)
  alpha<-c(d1, d2, 1-d1-d2)
  sigma1<-runif(K, 0.5, 0.8)
  sigma2<-runif(K, 0.5, 0.8)
  P<-matrix(rep(0,n*K),ncol=K)
  Q<-matrix(rep(0,n*K),ncol=K)
  R<-matrix(rep(0,n*K),ncol=K)
  H<-matrix(rep(0,n*K),ncol=K)
  ## Circle
  for (S in 1:200){
    ## E-step
    for (j in 1:K){
      for (i in 1:n){
        P[i,j]<-sapply(M[i], dnorm, beta1[j]+a[j]*X[i], sigma1[j])
        Q[i,j]<-sapply(Y[i], dnorm, beta2[j]+c[j]*X[i]+b[j]*M[i], sigma2[j])
        H[i,j]<-alpha[j]*P[i,j]*Q[i,j]
      }
    }
    H<-H/rowSums(H)
    oldbeta1<-beta1
    oldbeta2<-beta2
    olda<-a
    oldb<-b
    oldc<-c
    oldalpha<-alpha
    oldsigma1<-sigma1
    oldsigma2<-sigma2
    ## M-step
    for (j in 1:K){
      alpha[j]<-sum(H[,j])/sum(H)
      beta1[j]<-(solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M)[1]
      a[j]<-(solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M)[2]
      beta2[j]<-(solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y)[1]
      c[j]<-(solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y)[2]
      b[j]<-(solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y)[3]
      sigma1[j]<-{(M-c(beta1[j], a[j])%*%t(V))%*%diag(H[,j])%*%t(M-c(beta1[j], a[j])%*%t(V))/sum(H[,j])}^0.5
      sigma2[j]<-{(Y-c(beta2[j], c[j], b[j])%*%t(Z))%*%diag(H[,j])%*%t((Y-c(beta2[j], c[j], b[j])%*%t(Z)))/sum(H[,j])}^0.5
    }
    ## Change condition
    espsilo<-1e-5
    if (sum(abs(a-olda)<espsilo) &
        sum(abs(b-oldb)<espsilo) &
        sum(abs(c-oldc)<espsilo) &
        sum(abs(alpha-oldalpha)<espsilo)&
        sum(abs(beta1-oldbeta1)<espsilo)&
        sum(abs(beta2-oldbeta2)<espsilo)
    ) break
    cat('S', S, 'a', a, 'b', b, 'c', c,'alpha', alpha, '\n')
  }
  # Estimation
  men1<-which.min(c(a[1]*b[1], a[2]*b[2], a[3]*b[3]))
  men2<-which.max(c(a[1]*b[1], a[2]*b[2], a[3]*b[3]))
  men3<-setdiff(c(1,2,3), c(men1, men2))
  me1<-a[men1]*b[men1]
  me2<-a[men2]*b[men2]
  me3<-a[men3]*b[men3]
  
  a1_var<-sigma1[men1]^2*{solve(t(V)%*%diag(H[,men1])%*%V)%*%(t(V)%*%diag(H[,men1]^2)%*%V)%*%
      solve(t(V)%*%diag(H[,men1])%*%V)}[2,2]
  b1_var<-sigma2[men1]^2*{solve(t(Z)%*%diag(H[,men1])%*%Z)%*%(t(Z)%*%diag(H[,men1]^2)%*%Z)%*%
      solve(t(Z)%*%diag(H[,men1])%*%Z)}[3,3]
  a2_var<-sigma1[men2]^2*{solve(t(V)%*%diag(H[,men2])%*%V)%*%(t(V)%*%diag(H[,men2]^2)%*%V)%*%
      solve(t(V)%*%diag(H[,men2])%*%V)}[2,2]
  b2_var<-sigma2[men2]^2*{solve(t(Z)%*%diag(H[,men2])%*%Z)%*%(t(Z)%*%diag(H[,men2]^2)%*%Z)%*%
      solve(t(Z)%*%diag(H[,men2])%*%Z)}[3,3]
  a3_var<-sigma1[men3]^2*{solve(t(V)%*%diag(H[,men3])%*%V)%*%(t(V)%*%diag(H[,men3]^2)%*%V)%*%
      solve(t(V)%*%diag(H[,men3])%*%V)}[2,2]
  b3_var<-sigma2[men3]^2*{solve(t(Z)%*%diag(H[,men3])%*%Z)%*%(t(Z)%*%diag(H[,men3]^2)%*%Z)%*%
      solve(t(Z)%*%diag(H[,men3])%*%Z)}[3,3]
  
  ME3[f,]<-c(me1, me2, me3)
  A3[f,]<-c(a[men1], a[men2], a[men3])
  B3[f,]<-c(b[men1], b[men2], b[men3])
  ALPHA3[f,]<-c(alpha[men1], alpha[men2], alpha[men3])
}


## Mean and SD
XXX<-ME3[,2] ## XXX=A3[, 1:3], B3[, 1:3], ALPHA3[, 1:3], and ME3[, 1:3]. 
c(round(mean(XXX), digits=3), round(sd(XXX), digits=3))







