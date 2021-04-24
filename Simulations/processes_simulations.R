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
library(lpSolve)
library(glmnet)
library(MASS)

M = 100
K = 50

set.seed(3005)


##########################################################################
#### 1) VAR(1) sample paths simulations ##################################
##########################################################################

#' Import coefficients tables if on different session than the one
#' coefs_simulations.R' was run on
#' & restructure the tables as lists of lists
A.full = as.matrix(read.table('YOUR PATH/A.full.txt'))
B.full = as.matrix(read.table('YOUR PATH/B.full.txt'))

S = 4
A.list = vector(mode = 'list', length = S)
for(l in 1:S){
  A.list[[l]] = vector(mode = 'list', length = M)
  for(m in 1:M){
    A.list[[l]][[m]] = A.full[((l-1)*K*M + K*(m-1) + 1):((l-1)*K*M + K*m),]
  }
}

B.list = vector(mode = 'list', length = S)
for(l in 1:S){
  B.list[[l]] = B.full[,((l-1)*M+1):(l*M)]
}


#' 1 : LOW DIMENSIONAL  / NON-SPARSE
N = 200
X1.list = vector(mode = 'list', length = M)
sig = diag(1,K)

for(m in 1:M){
  X1.list[[m]] = VARMAsim(N, arlags = 1, phi = A.list[[1]][[m]], sigma = sig)$series
  colnames(X1.list[[m]]) = 1:K
}

#' 2 : HIGH DIMENSIONAL / NON-SPARSE
N = 40
X2.list = vector(mode = 'list', length = M)

for(m in 1:M){
  X2.list[[m]] = VARMAsim(N, arlags = 1, phi = A.list[[2]][[m]], sigma = sig)$series
  colnames(X2.list[[m]]) = 1:K
}

#' 3 : LOW DIMENSIONAL  / SPARSE
N = 200
X3.list = vector(mode = 'list', length = M)

for(m in 1:M){
  X3.list[[m]] = VARMAsim(N, arlags = 1, phi = A.list[[3]][[m]], sigma = sig)$series
  colnames(X3.list[[m]]) = 1:K
}

#' 4 : HIGH DIMENSIONAL / SPARSE
N = 40
X4.list = vector(mode = 'list', length = M)

for(m in 1:M){
  X4.list[[m]] = VARMAsim(N, arlags = 1, phi = A.list[[4]][[m]], sigma = sig)$series
  colnames(X4.list[[m]]) = 1:K
}


