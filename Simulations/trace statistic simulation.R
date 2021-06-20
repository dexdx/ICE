N = 2000                #number of observations
M = 5000                #number of iterations
D = 200                 #max number of system dimensions

jo.trace = matrix(0, nrow = M, ncol = D)

ptm = proc.time()


for(d in 1:D){
  
  for(m in 1:M){
    
    E = matrix(rnorm(N*d), ncol = d)
    X = apply(E, 2, cumsum)
    S = rbind(matrix(0, ncol = d), X)
    S = S[-(N+1),]
    
    sig.1 = t(E) %*% S
    sig.2 = t(S) %*% S
    sig.3 = t(S) %*% E
    
    mat = sig.1 %*% solve(sig.2) %*% sig.3
    mat.tr = sum(diag(mat))
    
    jo.trace[m,d] = mat.tr
    
  }
  
}


ptm = proc.time() - ptm
