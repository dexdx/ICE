ICE = function(PI, r){
  
  K = dim(PI)[1]
  PI.svd = svd(PI)
  PI_lr = PI.svd$u %*% diag(c(PI.svd$d[1:r],rep(0,K-r))) %*% t(PI.svd$v)
  PI_lr.rref = echelon(PI_lr)
  
  B = as.matrix(PI_lr.rref[1:r,])
  return(B)
  
}