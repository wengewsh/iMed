# Simulation with same direction mediation effects
library(latex2exp)
## NULL Setting
N<-1000
n<-500 # 100, 200, 300, 400, 500
ME<-matrix(0, N, 2)
A<-matrix(0, N, 2)
B<-matrix(0, N, 2)
ALPHA<-NULL

# Simulation 
for (f in 1:N)
{
  ## True parameters
  beta1<-c(1.0, -1.2)
  beta2<-c(0.6, -0.8)
  a<-c(0.4, 0.6)
  b<-c(0.5, 0.8)
  c<-c(-1.0, 1.8)
  sigma1<-c(1.2, 0.8)
  sigma2<-c(0.7, 1.3)
  alpha<-c(0.4,0.6)
  ## Sample 
  Y<-NULL
  M<-NULL
  n1<-floor(n*alpha[1])
  n2<-n-n1
  X<-rnorm(n)
  M[1:n1]<-beta1[1]+a[1]*X[1:n1]+rnorm(n1)*sigma1[1]
  Y[1:n1]<-beta2[1]+c[1]*X[1:n1]+b[1]*M[1:n1]+rnorm(n1)*sigma2[1]
  M[(1+n1):n]<-beta1[2]+a[2]*X[(1+n1):n]+rnorm(n2)*sigma1[2]
  Y[(1+n1):n]<-beta2[2]+c[2]*X[(1+n1):n]+b[2]*M[(1+n1):n]+rnorm(n2)*sigma2[2]
  U<-rep(1, n)
  V<-matrix(c(U, X), ncol=2)
  Z<-matrix(c(U, X, M), ncol=3)
  
  # EM for Two Mixed Gaussians
  ## Initialization parameters 
  K<-2
  beta1<-runif(K)
  beta2<-runif(K)
  a<-runif(K)
  b<-runif(K)
  c<-runif(K)
  d<-runif(1, 0.4, 0.6)
  alpha<-c(d,1-d)
  sigma1<-runif(K, 0.5, 1)
  sigma2<-runif(K, 0.5, 1)
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
  ## EM estimates
  men1<-which.min(c(a[1]*b[1],a[2]*b[2]))
  men2<-which.max(c(a[1]*b[1],a[2]*b[2]))
  me1<-a[men1]*b[men1]
  me2<-a[men2]*b[men2]
  
  a1_var<-sigma1[men1]^2*{solve(t(V)%*%diag(H[,men1])%*%V)%*%(t(V)%*%diag(H[,men1]^2)%*%V)%*%
      solve(t(V)%*%diag(H[,men1])%*%V)}[2,2]
  b1_var<-sigma2[men1]^2*{solve(t(Z)%*%diag(H[,men1])%*%Z)%*%(t(Z)%*%diag(H[,men1]^2)%*%Z)%*%
      solve(t(Z)%*%diag(H[,men1])%*%Z)}[3,3]
  a2_var<-sigma1[men2]^2*{solve(t(V)%*%diag(H[,men2])%*%V)%*%(t(V)%*%diag(H[,men2]^2)%*%V)%*%
      solve(t(V)%*%diag(H[,men2])%*%V)}[2,2]
  b2_var<-sigma2[men2]^2*{solve(t(Z)%*%diag(H[,men2])%*%Z)%*%(t(Z)%*%diag(H[,men2]^2)%*%Z)%*%
      solve(t(Z)%*%diag(H[,men2])%*%Z)}[3,3]
  
  ########## Interested Parameters' Estimation
  ME[f,]<-c(me1, me2)
  A[f,]<-c(a[men1], a[men2])
  B[f,]<-c(b[men1], b[men2])
  ALPHA[f]<-alpha[men1]
}


## Mean and SD of Interested Estimates 
XXX<-ME[,2]  ## XXX=A[, 1:2], B[, 1:2], ALPHA, and ME[, 1:2]. 
c(round(mean(XXX), digits=3), round(sd(XXX), digits=3))

