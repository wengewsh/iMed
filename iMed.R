#' @description The iMed function is the main function for detecting mediation effect subgroups and then estimate
 subgroup specific mediation effects.
#' @param A the exposure variable
#' @param X confounders
#' @param M the continuous mediator variable
#' @param Y the continuous outcome variable
#' @export

################################################################# iMed Function
iMed<-function(A, X, M, Y){

### NULL sBIC
sBIC<-NULL

### Needed parameter
L<-15
par<-length(X[1,])*2+7

#####################################################  K=1
### Observations
Data<-cbind(A, X, M, Y)  ## X is a matrix, not a vector.

# A is treatment
# X is a vector of covariates.
# M is a mediator
# Y is outcome

### New Data
n<-length(Y)
I<-rep(1, n)
V<-cbind(I, A, X)
Z<-cbind(I, M, A, X)


#### Estimation for K=1
a1<-(solve(t(V)%*%V)%*%t(V)%*%M)[2]
b1<-(solve(t(Z)%*%Z)%*%t(Z)%*%Y)[2]

## SD Estimation for K=1
sigma1_hat<-sqrt(sum((M-V%*%(solve(t(V)%*%V)%*%t(V)%*%M))^2)/n)
sigma2_hat<-sqrt(sum((Y-Z%*%(solve(t(Z)%*%Z)%*%t(Z)%*%Y))^2)/n)

## BICs for K=1
K<-1
P<-matrix(rep(0,n*K),ncol=K)
Q<-matrix(rep(0,n*K),ncol=K)
H<-matrix(rep(0,n*K),ncol=K)
for (i in 1:n){
  P[i]<-sapply(M[i], dnorm, V[i,]%*%(solve(t(V)%*%V)%*%t(V)%*%M), sigma1_hat)
  Q[i]<-sapply(Y[i], dnorm, Z[i,]%*%(solve(t(Z)%*%Z)%*%t(Z)%*%Y), sigma2_hat)
  H[i]<-P[i]*Q[i]
}

##### sBIC for K=1
L_11<-prod(H*L)/n^(par/2)
sBIC[1]<-log(L_11)-n*log(L)



############################################################ K=2
K<-2
num1<-length(X[1,])+2  ## It makes sense as X is a matrix.
num2<-num1+1
coefficients1<-matrix(rep(0,num1*K),ncol=K)
coefficients2<-matrix(rep(0,num2*K),ncol=K)
sigma1<-NULL
sigma2<-NULL
alpha2<-NULL
a2<-NULL
b2<-NULL

## Computation of Initialization parameters
H<-matrix(rep(0,n*K),ncol=K)
H[,1]<-runif(n, 0, 1)
H[,2]<-runif(n, 0, 1)

for (j in 1:K){
  alpha[j]<-sum(H[,j])/sum(H)
  coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
  coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
  sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M-V%*%coefficients1[,j])/sum(H[,j])}^0.5
  sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
}

### NULL matrix
P<-matrix(rep(0,n*K),ncol=K)
Q<-matrix(rep(0,n*K),ncol=K)
R<-matrix(rep(0,n*K),ncol=K)
H<-matrix(rep(0,n*K),ncol=K)

## Loops
for (S in 1:10000){
  ## E-step
  for (j in 1:K){
    for (i in 1:n){
      P[i, j]<-sapply(M[i], dnorm, V[i,]%*%coefficients1[,j], sigma1[j])
      Q[i, j]<-sapply(Y[i], dnorm, Z[i,]%*%coefficients2[,j], sigma2[j])
      R[i, j]<-alpha[j]*P[i,j]*Q[i,j]
    }
  }
  H<-R/rowSums(R)
  oldcoefficients1<-coefficients1
  oldcoefficients2<-coefficients2
  olda<-a
  oldb<-b
  oldalpha<-alpha
  oldsigma1<-sigma1
  oldsigma2<-sigma2
  ## M-step
  for (j in 1:K){
    alpha[j]<-sum(H[,j])/sum(H)
    coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
    coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
    a<-c(coefficients1[2, 1], coefficients1[2, 2])
    b<-c(coefficients2[2, 1], coefficients2[2, 2])
    sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M- V%*%coefficients1[,j])/sum(H[,j])}^0.5
    sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
  }
  ## Change condition
  espsilo<-1e-5
  if (sum(abs(a-olda)<espsilo) &
      sum(abs(b-oldb)<espsilo) &
      sum(abs(alpha-oldalpha)<espsilo)
  ) break
  cat('S', S, 'a', a, 'b', b, 'alpha', alpha, '\n')
}
### EM estimates
a2<-a
b2<-b
alpha2<-alpha

### sBIC for K=2
L_21<-prod((R[,1]+R[,2])*L)/n^((2*par)/2)
L_22<-prod((R[,1]+R[,2])*L)/n^((2*par+1)/2)
B<-L_11-L_22
C<-L_21*L_11
L_2<-(sqrt(B^2+4*C)-B)/2
sBIC[2]<-log(L_2)-n*log(L)


########################### K=3
K<-3
coefficients1<-matrix(rep(0,num1*K),ncol=K)
coefficients2<-matrix(rep(0,num2*K),ncol=K)
sigma1<-NULL
sigma2<-NULL
alpha3<-NULL
a3<-NULL
b3<-NULL

## Computation of Initialization parameters
H<-matrix(rep(0,n*K),ncol=K)
H[,1]<-runif(n, 0, 1)
H[,2]<-runif(n, 0, 1)
H[,3]<-runif(n, 0, 1)

for (j in 1:K){
  alpha[j]<-sum(H[,j])/sum(H)
  coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
  coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
  sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M-V%*%coefficients1[,j])/sum(H[,j])}^0.5
  sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
}

### NULL matrix
P<-matrix(rep(0,n*K),ncol=K)
Q<-matrix(rep(0,n*K),ncol=K)
R<-matrix(rep(0,n*K),ncol=K)
H<-matrix(rep(0,n*K),ncol=K)

## Loops
for (S in 1:10000){
  ## E-step
  for (j in 1:K){
    for (i in 1:n){
      P[i, j]<-sapply(M[i], dnorm, V[i,]%*%coefficients1[,j], sigma1[j])
      Q[i, j]<-sapply(Y[i], dnorm, Z[i,]%*%coefficients2[,j], sigma2[j])
      R[i, j]<-alpha[j]*P[i,j]*Q[i,j]
    }
  }
  H<-R/rowSums(R)
  oldcoefficients1<-coefficients1
  oldcoefficients2<-coefficients2
  olda<-a
  oldb<-b
  oldalpha<-alpha
  oldsigma1<-sigma1
  oldsigma2<-sigma2
  ## M-step
  for (j in 1:K){
    alpha[j]<-sum(H[,j])/sum(H)
    coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
    coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
    a<-c(coefficients1[2, 1], coefficients1[2, 2])
    b<-c(coefficients2[2, 1], coefficients2[2, 2])
    sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M- V%*%coefficients1[,j])/sum(H[,j])}^0.5
    sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
  }
  ## Change condition
  espsilo<-1e-5
  if (sum(abs(a-olda)<espsilo) &
      sum(abs(b-oldb)<espsilo) &
      sum(abs(alpha-oldalpha)<espsilo)
  ) break
  cat('S', S, 'a', a, 'b', b, 'alpha', alpha, '\n')
}
### EM estimates
a3<-a
b3<-b
alpha3<-alpha

### sBIC for K=3
L_31<-prod((R[,1]+R[,2]+R[,3])*L)/n^((3*par+0)/2)
L_32<-prod((R[,1]+R[,2]+R[,3])*L)/n^((3*par+1)/2)
L_33<-prod((R[,1]+R[,2]+R[,3])*L)/n^((3*par+2)/2)
B<-L_11+L_2-L_33
C<-L_31*L_11+L_32*L_2
L_3<-(sqrt(B^2+4*C)-B)/2
sBIC[3]<-log(L_3)-n*log(L)

########################### K=4
K<-4
coefficients1<-matrix(rep(0,num1*K),ncol=K)
coefficients2<-matrix(rep(0,num2*K),ncol=K)
sigma1<-NULL
sigma2<-NULL
alpha4<-NULL
a4<-NULL
b4<-NULL

## Computation of Initialization parameters
H<-matrix(rep(0,n*K),ncol=K)
H[,1]<-runif(n, 0, 1)
H[,2]<-runif(n, 0, 1)
H[,3]<-runif(n, 0, 1)
H[,4]<-runif(n, 0, 1)

for (j in 1:K){
  alpha[j]<-sum(H[,j])/sum(H)
  coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
  coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
  sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M-V%*%coefficients1[,j])/sum(H[,j])}^0.5
  sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
}

### NULL matrix
P<-matrix(rep(0,n*K),ncol=K)
Q<-matrix(rep(0,n*K),ncol=K)
R<-matrix(rep(0,n*K),ncol=K)
H<-matrix(rep(0,n*K),ncol=K)

## Loops
for (S in 1:10000){
  ## E-step
  for (j in 1:K){
    for (i in 1:n){
      P[i, j]<-sapply(M[i], dnorm, V[i,]%*%coefficients1[,j], sigma1[j])
      Q[i, j]<-sapply(Y[i], dnorm, Z[i,]%*%coefficients2[,j], sigma2[j])
      R[i, j]<-alpha[j]*P[i,j]*Q[i,j]
    }
  }
  H<-R/rowSums(R)
  oldcoefficients1<-coefficients1
  oldcoefficients2<-coefficients2
  olda<-a
  oldb<-b
  oldalpha<-alpha
  oldsigma1<-sigma1
  oldsigma2<-sigma2
  ## M-step
  for (j in 1:K){
    alpha[j]<-sum(H[,j])/sum(H)
    coefficients1[, j]<-solve(t(V)%*%diag(H[,j])%*%V)%*%t(V)%*%diag(H[,j])%*%M
    coefficients2[, j]<-solve(t(Z)%*%diag(H[,j])%*%Z)%*%t(Z)%*%diag(H[,j])%*%Y
    a<-c(coefficients1[2, 1], coefficients1[2, 2])
    b<-c(coefficients2[2, 1], coefficients2[2, 2])
    sigma1[j]<-{t(M-V%*%coefficients1[,j])%*%diag(H[,j])%*%(M- V%*%coefficients1[,j])/sum(H[,j])}^0.5
    sigma2[j]<-{t(Y-Z%*%coefficients2[,j])%*%diag(H[,j])%*%(Y-Z%*%coefficients2[,j])/sum(H[,j])}^0.5
  }
  ## Change condition
  espsilo<-1e-5
  if (sum(abs(a-olda)<espsilo) &
      sum(abs(b-oldb)<espsilo) &
      sum(abs(alpha-oldalpha)<espsilo)
  ) break
  cat('S', S, 'a', a, 'b', b, 'alpha', alpha, '\n')
}

### EM estimates
a4<-a
b4<-b
alpha4<-alpha

### sBIC for K=4
L_41<-prod((R[,1]+R[,2]+R[,3]+R[,4])*L)/n^((4*par+0)/2)
L_42<-prod((R[,1]+R[,2]+R[,3]+R[,4])*L)/n^((4*par+1)/2)
L_43<-prod((R[,1]+R[,2]+R[,3]+R[,4])*L)/n^((4*par+2)/2)
L_44<-prod((R[,1]+R[,2]+R[,3]+R[,4])*L)/n^((4*par+3)/2)
B<-L_11+L_2+L_3-L_44
C<-L_41*L_11+L_42*L_2+L_43*L_3
L_4<-(sqrt(B^2+4*C)-B)/2
sBIC[4]<-log(L_4)-n*log(L)

K_opt<-which.max(sBIC)

Parameters<-NULL

if (length(a1)==K_opt){
  ME<-a1*b1
  Parameters<-list(a1, b1, ME, 1)
} else if (length(a2)==K_opt){
  ME<-a2*b2
  Parameters<-list(a2, b2, ME, alpha2)
  } else if (length(a3)==K_opt){
    ME<-a3*b3
    Parameters<-list(a3, b3, ME, alpha3)
  } else {
    ME<-a4*b4
    Parameters<-list(a4, b4, ME, alpha4)
  }
results<-list(K_opt=K_opt, a=Parameters[1], b=Parameters[2], ME=Parameters[3], alpha=Parameters[4])
return(results)
}


