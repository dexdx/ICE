##########################################
##### FONCTIONS UTILITAIRES DIVERSES #####
##########################################

## Returns array index of lower triangle elements excluding diagonal
lowtri.idx = function(K){
  for(i in 1:(K-1)){
    if(i == 1){
      idx = (i*K - (K-i-1)):(i*K)
    } else{
      idx = c(idx, (i*K - (K-i-1)):(i*K))
    }
  }
  return(idx)
}

## Translates array index into matrix index (row,col) for square matrix
ArrToMat.idx = function(arr.idx, K){
  mat.idx = numeric(2)
  if(arr.idx %% K == 0){
    mat.idx = c(K, arr.idx/K + 1)
  } else{
    j = floor(arr.idx/K) + 1
    if(j == 1){
      i = arr.idx
    } else{
      i = arr.idx - K*(j-1)
    }
    mat.idx = c(i,j)
  }
  return(mat.idx)
}

## Translates matrix index into array index for square matrix
MatToArr.idx = function(mat.idx, K){
  return( (mat.idx[2]-1)*K + mat.idx[1] )
}

## Inverse vector operator
VecToMat = function(y,n.row,n.col){
  y = as.numeric(y)
  n.vec = length(y)
  if(n.vec == prod(n.row,n.col)){
    M = matrix(0,n.row,n.col)
    for(j in 1:n.col){
      M[,j] = y[((j-1)*n.row + 1):(j*n.row)]
    }
    return(M)
  } else{
    print('Error: incorrect dimensions provided')
    return(NULL)
  }
}

## Returns n fitted values of VAR(1) for initial values and coefficient matrix
VAR1.fit = function(initial, n, A){
  x.fit = matrix(0,n,ncol(A))
  initial = t(as.matrix(as.numeric(initial)))
  x.fit[1,] = initial
  for(t in 2:n){
    if(t == 2){
      x.fit[t,] = initial %*% t(A)
    } else{
      x.fit[t,] = t(as.matrix(x.fit[t-1,])) %*% t(A)
    }
  }
  return(x.fit)
}

## Not in operator
`%notin%` <- Negate(`%in%`)

## Tries to compute true inverse, if non feasible computes pseudo inverse
forceInv = function(x){
  error = class(try(solve(x), silent = T))
  if("try-error" %notin% error){
    return(solve(x))
  } else return(ginv(x))
}

###################################
##### FONCTIONS D'ESTIMATEURS #####
###################################

##### Multivariate Lasso Estimator
MVlasso = function(x){
  
  N = dim(x)[1]
  K = dim(x)[2]
  Y = x[-1,]
  X = x[-N,]
  
  net = glmnet(X, Y, family = 'mgaussian')
  num = length(net$lambda)
  
  A.net.list = vector(mode = 'list', length = num)
  for(m in 1:num){
    for(k in 1:K){
      if(k == 1){
        A.net.list[[m]] = t(as.matrix(net$beta[[k]][,m]))
      } else{
        A.net.list[[m]] = rbind(A.net.list[[m]], t(as.matrix(net$beta[[k]][,m])))
      }
    }
    for(k in 1:K){
      if(sum(A.net.list[[m]][,k]==0)==K) A.net.list[[m]][k,k] = 1
    }
  }

  return(A.net.list)
  
}


##### Sparse Yule-Waler Estimator
sYWE = function(x){
  
  N = dim(x)[1]
  K = dim(x)[2]
  
  S = crossprod(x)/N
  
  X0 = x[-1,]
  X1 = x[-N,]
  
  S1 = crossprod(X0,X1)/(N-1)
  
  # lam0 = max(abs((matrix(rnorm(K^2),K) %*% S) - S1))
  # lams = seq(0, 2*lam0, (2*lam0)/(50-1))
  
  lams = seq(0,1,1/50)
  l = length(lams)
  
  lp.lam.list = vector(mode = 'list', length = l)
  lp.k.list = vector(mode = 'list', length = K)
  
  for(lam in 1:l){
    
    A.lam = matrix(0,K,K)
    
    for(k in 1:K){
      
      f.obj = rep(1, 2*K)
      f.con1 = cbind(-S,S)
      f.con2 = cbind(S,-S)
      f.con3 = diag(1,2*K)
      f.con = rbind(f.con1,f.con2,f.con3)
      f.dir = rep('>=', 4*K)
      f.rhs = c( -S1[k,] - lams[lam], S1[k,] - lams[lam], rep(0,2*K) )
      
      lp.k.list[[k]] = lp('min', f.obj, f.con, f.dir, f.rhs)
      A.k = lp.k.list[[k]][["solution"]][1:K] - lp.k.list[[k]][["solution"]][(K+1):(2*K)]
      A.lam[k,] = A.k
      
    }
    
    for(k in 1:K){
      if(sum(A.lam[,k]==0)==K) A.lam[k,k] = 1
    }
    
    lp.lam.list[[lam]] = A.lam
    
  }

  return(lp.lam.list)
  
}


## Estimates PSC matrix for a given frequency (w) and kernel coefficient (v)
PSC = function(x, w, v, n.min=10){
  
  N = dim(x)[1]
  K = dim(x)[2]
  
  l.max = (N - n.min)/2
  im = complex(1,0,1)
  
  Covs = vector(mode = 'list', length = l.max+1)
  
  for(l in 0:l.max){
    if(l == 0){
      Covs[[l+1]] = crossprod(x)
    } else{
      X1 = x[-(1:l),]
      X2 = x[-((N-l+1):N),]
      Covs[[l+1]] = crossprod(X1,X2)
    }
  }
  
  for(l in 1:l.max){
    
    kern = sin(v*l)/(v*l)
    
    if(l == 1){
      SiCov = kern * ( Covs[[l+1]] * exp(-im*w*l) + t(Covs[[l+1]]) *exp(-im*w*l) )
    } else{
      SiCov = SiCov + kern * ( Covs[[l+1]] * exp(-im*w*l) + t(Covs[[l+1]]) *exp(-im*w*l) )
    }
    
  }
  
  Sx = 1/(2*pi) * (Covs[[1]] + SiCov)
  
  return(Mod(solve(Sx)))
  
}


###### Sparse PSC Estimator
sPSC = function(x){
  
  N = dim(x)[1]
  K = dim(x)[2]
  
  freqs = numeric(N)
  for(i in 1:N) freqs[i] = (2*pi*i)/N
  
  pair_ind = lowtri.idx(K)
  pair.len = length(pair_ind)
  freq.len = length(freqs)
  
  pairPSC.w = matrix(0,pair.len,freq.len)
  
  for(w in 1:freq.len){
    
    PSC.w = PSC(x,freqs[w],1)
    pairPSC.w[,w] = PSC.w[pair_ind]
    
  }
  
  PSC.max = apply(pairPSC.w, 1, max)
  PSC.max.ord = pair_ind[order(PSC.max)]
  
  
  Ms = seq(0,pair.len,floor(pair.len/50))
  Ms = Ms[-1]
  Ms.l = length(Ms)
  A.m.list = vector(mode = 'list', length = Ms.l)
  
  for(m in 1:Ms.l){
    
    sparse.idx = rep(1, K^2)
    lowtri.sparse = PSC.max.ord[1:Ms[m]]
    sparse.idx[lowtri.sparse] = 0
    
    lowtri.sparse.mat = sapply(lowtri.sparse, function(x) ArrToMat.idx(x,K))
    uptri.sparse.mat = lowtri.sparse.mat[c(2,1),]
    uptri.sparse = apply(uptri.sparse.mat, 2, function(x) MatToArr.idx(x,K))
    sparse.idx[uptri.sparse] = 0
    
    nz = sum(sparse.idx==1)
    nz.idx = which(sparse.idx == 1)
    R = matrix(0,K^2,nz)
    for(i in 1:nz){
      R[nz.idx[i],i] = 1
    }
    
    Y = x[-1,]
    X = x[-N,]
    
    S0 = crossprod(x)/N
    S0.inv = forceInv(S0)
    KR1 = kronecker(crossprod(X),S0.inv)
    KR2 = kronecker(t(X),S0.inv)
    cross1 = crossprod(R, KR1)
    cross2 = crossprod(R, KR2)
    inv2 = cross1 %*% R
    inv2 = forceInv(inv2)
    M1 = R %*% inv2
    M2 = cross2 %*% as.matrix(as.numeric(t(Y)))
    phi = M1 %*% M2
    A = VecToMat(phi,K,K)
    
    X.fit = VAR1.fit(X[1,], nrow(X), A)
    res.fit = X - X.fit
    
    S1 = crossprod(res.fit)/(N-1)
    S.dif0 = sum(abs(S1 - S0))
    S.dif1 = S.dif0 - 1
    dif.tol = 100
    tours = 0
    
    error = class(try(
    while(S.dif1 > dif.tol && S.dif1 < S.dif0){
      tours = tours + 1
      S0 = S1
      S.dif0 = S.dif1
      
      S0.inv = forceInv(S0)
      KR1 = kronecker(crossprod(X),S0.inv)
      KR2 = kronecker(t(X),S0.inv)
      cross1 = crossprod(R, KR1)
      cross2 = crossprod(R, KR2)
      inv2 = cross1 %*% R
      inv2 = forceInv(inv2)
      M1 = R %*% inv2
      M2 = cross2 %*% as.matrix(as.numeric(t(Y)))
      phi = M1 %*% M2
      
      A = VecToMat(phi,K,K)
      X.fit = VAR1.fit(X[1,], nrow(X), A)
      res.fit = X - X.fit
      S1 = crossprod(res.fit)/(N-1)
      S.dif1 = sum(abs(S1 - S0))
      
    }
    , silent = T))
    
    if(error == "try-error"){
      A.m.list[[m]] = diag(1,K)
    } else{
      A.m.list[[m]] = A
    }
      
  
  }
  
  return(A.m.list)
  
}




