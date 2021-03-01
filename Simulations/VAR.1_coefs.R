#' ALL VARs ARE "NEARLY" INTEGRATED (non-trivial cointegrating relationships)

#' 1 : LOW DIMENSIONAL  / NON-SPARSE
#' 2 : HIGH DIMENSIONAL / NON-SPARSE
#' 3 : LOW DIMENSIONAL  / SPARSE
#' 4 : HIGH DIMENSIONAL / SPARSE

#' Both low and high dimensional simulated processes are set up with K = 50 as the number
#' of variables. Only the number of observations changes between these two settings, with
#' the low dimensional setting being comprised of N = 40 observations and the high dimensional
#' one of N = 200.

#' Each of the 4 set ups is simulated M = 100 times.


library(MTS)
library(urca)
library(Matrix)
library(matlib)

M = 100

set.seed(3005)


##########################################################################
#### 1) Non-Explosive Transition Matrices Simulation #####################
##########################################################################

K = 50
S = 4
A.list = vector(mode = 'list', length = S)
for(l in 1:S) A.list[[s]] = vector(mode = 'list', length = M)

sig = diag(1,K)
diag.coefs = c(0.8,1)

for(l in 1:S){
  for(m in 1:M){
    
    A = matrix(0,K,K)
    expl = T
    
    while(expl == T){
      for(k in 1:K){
        A[k,k] = runif(1,diag.coefs[1],diag.coefs[2])
        reste = 1 - A[k,k]
        A[k,-k] = runif(K-1, -2*reste/(K-1), 2*reste/(K-1))
      }
      eigs = Mod(eigen(A)$values)
      if(max(eigs) >= 1){
        expl = T
      } else expl = F
    }
    
    A.list[[l]][[m]] = A
    
  }
}


##########################################################################
#### 2) - Transformation into low rank (r=1) approximated matrices #######
####    - Retrieve resulting true cointegrating vector             #######
##########################################################################

B.list = vector(mode = 'list', length = S)
for(l in 1:S) B.list[[s]] = vector(mode = 'list', length = M)

## Non-Sparse cases (1 and 2)
for(l in 1:2){
  for(m in 1:M){
    
    A = A.list[[l]][[m]]
    PI = A - diag(1,K)
    
    # PI low rank approximation
    PI.svd = svd(PI)
    r = 1
    PI_lr = PI.svd$u %*% diag(c(PI.svd$d[1:r],rep(0,K-r))) %*% t(PI.svd$v)
    A.lr = PI_lr + diag(1,K)
    
    # PI full rank factorization
    PI_lr.rref = echelon(PI_lr)
    C = as.matrix(PI_lr[,1:r])
    B = t(PI_lr.rref[1:r,])
    
    # Modify list
    A.list[[l]][[m]] = A.lr
    
    # Save cointegrating vector
    B.list[[l]][[m]] = B
    
  }
}

## Sparse cases (3 and 4)
for(l in 3:4){
  for(m in 1:M){
    
    A = A.list[[l]][[m]]
    PI = A - diag(1,K)
    
    # PI low rank approximation
    PI.svd = svd(PI)
    r = 1
    PI_lr = PI.svd$u %*% diag(c(PI.svd$d[1:r],rep(0,K-r))) %*% t(PI.svd$v)
    A.lr = PI_lr + diag(1,K)
    
    # PI full rank factorization
    PI_lr.rref = echelon(PI_lr)
    C = as.matrix(PI_lr[,1:r])
    B = t(PI_lr.rref[1:r,])
    
    # Induce (i.i.d.) sparsity in cointegrating vector
    B.s = matrix(0,nrow=1,ncol=K)
    nz = which(rbinom(K,1,1/2)==1)
    B.s[nz] = B[nz]
    
    # Rewrite sparse transition matrix
    PI_lr.s = C %*% B.s
    A.s = PI_lr.s + diag(1,K)
    
    # Modify list
    A.list[[l]][[m]] = A.s
    
    # Save (sparse) cointegrating vector
    B.list[[l]][[m]] = B.s
    
  }
}


##########################################################################
#### 3) Save lists in tables for use in process simulations file #########
##########################################################################

#' This step is only necessary if the 'process simulations' file is run on another
#' session in which case the above generated lists lie not in the global environment.

# Stack matrices of lists for saving to table file
for(m in 1:M){
  if(m == 1){
    A.full = A.list[[1]]
    B.full = t(B.list[[1]])
  } else{
    A.full = rbind(A.full,A.list[[m]])
    B.full = cbind(B.full, t(B.list[[m]]))
  }
}






