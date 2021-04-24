#########################################################################################
# FUNCTIONS TO IMPLEMENT THE SPARSE COINTEGRATION METHOD AS DESCRIBED IN                #
# Wilms I. and CROUX C. (2016), "Forecasting using sparse cointegration." Working paper #
#########################################################################################

# REQUIRED LIBRARIES
library(mvtnorm)
library(glmnet)
library(MASS)
library(glasso)
library(lars)
library(huge)


SparseCointegration_Lasso<-function(p,Y,X,Z,r,alpha.init=NULL,beta.init=NULL,max.iter=10,conv=10^-2,lambda.gamma=matrix(seq(from=0.01,to=0.001,length=5),nrow=1),rho.glasso=seq(from=1,to=0.1,length=5),lambda_beta=matrix(seq(from=0.1,to=0.001,length=100),nrow=1),cutoff=0.8,glmnetthresh=1e-04){
  ### Main function to perform sparse cointegration ###
  
  ## INPUT
  # p: number of lagged differences
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r: cointegration rank
  # alpha.init: initial value for adjustment coefficients
  # beta.init: initial value for cointegrating vector
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # lambda.gamma: tuning paramter short-run effects
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # glmnetthresh: tolerance parameter glmnet function
  
  ## OUTPUT
  # BETAhat: estimate of cointegrating vectors 
  # ALPHAhat: estimate of adjustment coefficients
  # ZBETA: estimate of short-run effects
  # OMEGA: estimate of inverse covariance matrix

  
  ## START CODE
  
  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  # Starting value
  if(is.null(alpha.init) & is.null(beta.init)){
    SigmaYY<-diag(apply(Y,2,var))
    SigmaZZ<-diag(apply(Z,2,var))
    
    SigmaYZ<-cov(Y,Z)
    SigmaZY<-cov(Z,Y)
    
    matrix<-solve(SigmaZZ)%*%SigmaZY%*%solve(SigmaYY)%*%SigmaYZ
    decomp<-eigen(matrix)
    beta.init<-decomp$vectors[,1:r]
    
    SVDinit<-svd(t(Z%*%beta.init)%*%(Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T)))
    alpha.init<-t(SVDinit$u%*%t(SVDinit$v))
    beta.init<-matrix(NA,ncol=r,nrow=q)
    Yinit<-Re((Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T))%*%alpha.init)
    for (i in 1:r){
      FIT<-lars(x=Z,y=Yinit[,i],type="lasso",normalize=F,intercept=F)$beta[-1,]
      beta.init[,i]<-FIT[nrow(FIT),]
    } 
    Pi.init<-alpha.init%*%t(beta.init)
  }
  
  if(is.null(alpha.init) & !is.null(beta.init)){
    SVDinit<-svd(t(Z%*%beta.init)%*%(Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T)))
    alpha.init<-t(SVDinit$u%*%t(SVDinit$v))
    beta.init<-matrix(NA,ncol=r,nrow=q)
    Yinit<-(Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T))%*%alpha.init
    for (i in 1:r){
      FIT<-lars(x=Z,y=Yinit[,i],type="lasso",normalize=F,intercept=F)$beta[-1,]
      beta.init[,i]<-FIT[nrow(FIT),]
    } 
    Pi.init<-alpha.init%*%t(beta.init)
  }
  
  if(!is.null(alpha.init) & !is.null(beta.init)){
    Pi.init<-alpha.init%*%t(beta.init)
  }
  
  # Convergence parameters: initialization
  it<-1
  diff.obj<-10*conv
  Omega.init=diag(1,q)
  value.obj<-matrix(NA,ncol=1,nrow=max.iter+1)
  RESID<-Y-X%*%matrix(rep(diag(1,q),p-1),ncol=q,byrow=T)-Z%*%beta.init%*%t(alpha.init)
  value.obj[1,]<-(1/n)*sum(diag(RESID%*%Omega.init%*%t(RESID)))-log(det(Omega.init))
  
  
  while( (it<max.iter) &  (diff.obj>conv) )
  {
    # Obtain Gamma
    FIT1<-NTS.GAMMA.LassoReg(Y=Y,X=X,Z=Z,Pi=Pi.init,p=p,Omega=Omega.init,lambda.gamma=lambda.gamma,glmnetthresh=glmnetthresh)
    # Obtain Alpha
    FIT2<-NTS.ALPHA.Procrusted(Y=Y,X=X,Z=Z,ZBETA=FIT1$ZBETA,r=r,Omega=Omega.init,P=FIT1$P,beta=beta.init)
    # Obtain Beta and Omega
    FIT3<-NTS.BETA(Y=Y,X=X,Z=Z,ZBETA=FIT1$ZBETA,r=r,Omega=Omega.init,P=FIT1$P,alpha=FIT2$ALPHA,alphastar=FIT2$ALPHAstar,lambda_beta=lambda_beta,rho.glasso=rho.glasso,cutoff=cutoff,glmnetthresh=glmnetthresh)

    
    # Check convergence
    beta.new<-FIT3$BETA
    RESID<-Y-X%*%FIT1$ZBETA-Z%*%FIT3$BETA%*%t(FIT2$ALPHA)
    value.obj[1+it,]<- (1/n)* sum(diag((RESID)%*%FIT3$OMEGA%*%t(RESID))) - log(det(FIT3$OMEGA))
    diff.obj<-abs(value.obj[1+it,]-value.obj[it,])/abs(value.obj[it,])
    beta.init<-matrix(beta.init,nrow=q,ncol=r)
    alpha.init<-FIT2$ALPHA
    beta.init<-beta.new
    Pi.init<-alpha.init%*%t(beta.new) 
    Omega.init<-FIT3$OMEGA
    it<-it+1
  }
  
  out<-list(BETAhat=FIT3$BETA,ALPHAhat=FIT2$ALPHA,it=it,ZBETA=FIT1$ZBETA,OMEGA=FIT3$OMEGA)
}

NTS.GAMMA.LassoReg<-function(Y,X,Z,Pi,Omega,lambda.gamma=matrix(seq(from=1,to=0.01,length=10),nrow=1),p,cutoff=0.8,intercept=F,glmnetthresh=1e-04){
  ### FUNCTIONS TO ESTIMATE GAMMA, using L1 PENALTY ###
  
  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # Pi: estimate of cointegration space
  # Omega: estimate of inverse error covariance matrix
  # lambda.gamma: tuning paramter short-run effects
  # p: number of lagged differences
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept
  # glmnetthresh: tolerance parameter glmnet function
  
  ## OUTPUT
  # ZBETA: estimate of short-run effects
  # P: transformation matrix P derived from Omega
  
  
  ## START CODE
  
  # Setting dimensions
  q<-ncol(Y)
  
  # Creating data matrices
  Y.matrix<-Y-Z%*%t(Pi)
  X.matrix<-cbind(1,X)
  if(intercept==T){
    X.matrix<-cbind(1,X)
  } else {
    X.matrix<-cbind(X)
  }
  
  # Compute transformation matrix P
  decomp<-eigen(Omega)
  P<-(decomp$vectors)%*%diag(sqrt(decomp$values))%*%solve(decomp$vectors)
  
  
  # Final estimate 
  Pfinal<-kronecker(P,diag(rep(1,nrow(X.matrix))))   
  YmatrixP<-Pfinal%*%c(Y.matrix)
  YmatrixP.sd<-sd(Y.matrix)
  YmatrixP.scaled<-YmatrixP/YmatrixP.sd
  XmatrixP<-kronecker(P%*%diag(1,q),diag(rep(1,nrow(X.matrix)))%*%X.matrix)
  
  # Lasso Estimation
  # Determine lambda sequence: exclude all zero-solution
  determine_lambdasequence<-glmnet(y=YmatrixP.scaled,x=XmatrixP,standardize=F,intercept=F,lambda=lambda.gamma,family="gaussian",thresh=glmnetthresh)
  lambda_restricted<-matrix(lambda.gamma[1,which(determine_lambdasequence$df!=0)],nrow=1)
  
  if(length(which(determine_lambdasequence$df!=0))<=1){
    determine_lambdasequence<-glmnet(y=YmatrixP.scaled,x=XmatrixP,standardize=F,intercept=F,family="gaussian",thresh=glmnetthresh)
    lambda_restricted<-matrix(determine_lambdasequence$lambda,nrow=1)
  }
  
  # Time series cross-validation to determine value of the tuning parameter lambda_beta
  cutoff.n<-round(cutoff*nrow(Y))
  CVscore.gamma<-matrix(NA,ncol=ncol(lambda_restricted),nrow=nrow(Y)-cutoff.n)

  for (i in cutoff.n:nrow(Y)-1){ #Loop to calculate cross-validation score
    # Training Data
    Ytrain<-YmatrixP[1:i,]
    Ytrain.sd<-sd(Ytrain)
    Ytrain.scaled<-Ytrain/Ytrain.sd
    Xtrain<-XmatrixP[1:i,]
    
    # Test Data
    Ytest<-YmatrixP[i+1,]
    Xtest<-XmatrixP[i+1,]
    
    # Estimation
    GAMMA.scaled<-glmnet(y=Ytrain.scaled,x=Xtrain,lambda=lambda_restricted,standardize=F,intercept=F,family="gaussian",thresh=glmnetthresh)
    G_GAMMASD<-matrix(GAMMA.scaled$beta[,which(GAMMA.scaled$df!=0)],nrow=ncol(Xtrain)) # GAMMA in standardized scale
    G_GAMMA<-as.matrix(G_GAMMASD*Ytrain.sd)
    
    G_CV<-apply(G_GAMMA,2,CVscore.Reg,Y.data=Ytest,X.data=Xtest,Y.sd=Ytrain.sd) # CVscore
    CVscore.gamma[i-cutoff.n+1,c(which(GAMMA.scaled$df!=0))]<-G_CV  
  }

  CVscore.AVAILABLE<-as.matrix(CVscore.gamma[,apply(CVscore.gamma,2,AVAILABLE_LAMBDA)])
  lambda_restricted_AVAILABLE<-lambda_restricted[apply(CVscore.gamma,2,AVAILABLE_LAMBDA)]
  lambda.opt<-lambda_restricted_AVAILABLE[which.min(apply(CVscore.AVAILABLE,2,AVERAGE))]
  

  LASSOfinal<-glmnet(y=YmatrixP.scaled,x=XmatrixP,standardize=F,intercept=F,lambda=lambda.opt,family="gaussian",thresh=glmnetthresh)
  GAMMAscaled<-matrix(LASSOfinal$beta,ncol=1)
  GAMMA.Sparse<-GAMMAscaled*YmatrixP.sd
  
  
  ZBETA<-matrix(NA,ncol=q,nrow=ncol(X.matrix))
  for (i.q in 1:q)
  {
    if(intercept==T){
      ZBETA[,i.q]<-GAMMA.Sparse[((i.q-1)*((p-1)*q+1)+1):(i.q*(1+(p-1)*q))]
    } else {
      ZBETA[,i.q]<-GAMMA.Sparse[((i.q-1)*((p-1)*q)+1):(i.q*((p-1)*q))]
    }  
  }  
  
  out<-list(ZBETA=ZBETA,P=P)
}

NTS.ALPHA.Procrusted<-function(Y,X,Z,ZBETA,r,Omega,P,beta,intercept=F){
  ### FUNCTION TO ESTIMATE ALPHA ###
  
  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels 
  # ZBETA: estimate of short-run effects
  # r: cointegration rank
  # Omega: estimate of inverse error covariance matrix
  # P: transformation matrix P derived from Omega
  # beta: estimate of cointegrating vector
  # intercept: F do not include intercept, T include intercept in estimation short-run effects
  
  ## OUTPUT
  # ALPHA: estimate of adjustment coefficients
  # ALPHAstar: estimate of transformed adjustment coefficients
  
  ## START CODE 
  # Data matrices
  if(intercept==T){
    Ymatrix<-Y-cbind(1,X)%*%ZBETA
  } else {
    Ymatrix<-Y-cbind(X)%*%ZBETA
  }
  
  Xmatrix<-Z%*%beta
  decomp<-eigen(Omega)
  Pmin<-(decomp$vectors)%*%diag(1/sqrt(decomp$values))%*%solve(decomp$vectors) 
  
  # Singular Value Decomposition to compute ALPHA
  SingDecomp<-svd(t(Xmatrix)%*%Ymatrix%*%P)
  ALPHAstar<-t(SingDecomp$u%*%t(SingDecomp$v))
  if (r==1){
    ALPHAstar<-matrix(ALPHAstar,ncol=1)
  }
  ALPHA<-Pmin%*%ALPHAstar

  out<-list(ALPHA=ALPHA,ALPHAstar=ALPHAstar)
}

NTS.BETA<-function(Y,X,Z,ZBETA,r,Omega,P,alpha,alphastar,lambda_beta=matrix(seq(from=2,to=0.001,length=100),nrow=1),rho.glasso,cutoff,intercept=F,glmnetthresh=1e-04){
  ### FUNCTIONS TO ESTIMATE BETA AND OMEGA ###
  # First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
  # Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso 
  
  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels 
  # ZBETA: estimate of short-run effects
  # r: cointegration rank
  # Omega: estimate of inverse error covariance matrix
  # P: transformation matrix P derived from Omega
  # alpha: estimate of adjustment coefficients
  # alphastar: estimate of transformed adjustment coefficients
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept in estimation short-run effects
  # glmnetthresh: tolerance parameter glmnet function
  
  ## OUTPUT
  # BETA: estimate of cointegrating vectors 
  # OMEGA: estimate of inverse covariance matrix
  
  ## START CODE
  
  # Data matrices
  if(intercept==T){
    Ymatrix<-(Y-cbind(1,X)%*%ZBETA)%*%t(P)%*%alphastar
  } else {
    Ymatrix<-(Y-cbind(X)%*%ZBETA)%*%t(P)%*%alphastar
  }
  
  Xmatrix<-Z
  n<-nrow(Ymatrix)
  
  # Store Estimates of cointegrating vector
  BETA.Sparse<-matrix(NA,ncol=r,nrow=ncol(Y))
  
  
  # Performe Lasso
  for (i.r in 1:r)
  { # Determine each cointegrating vector by a Lasso Regression
    
    # Standardized Response
    Ymatrixsd<-Ymatrix[,i.r]/sd(Ymatrix[,i.r])
    
    # Determine lambda sequence: exclude all zero-solution
    determine_lambdasequence<-glmnet(y=Ymatrixsd,x=Xmatrix,standardize=F,intercept=F,lambda=lambda_beta,family="gaussian",thresh=glmnetthresh)
    lambda_restricted<-matrix(lambda_beta[1,which(determine_lambdasequence$df!=0)],nrow=1)
    if(length(which(determine_lambdasequence$df!=0))<=1){
      determine_lambdasequence<-glmnet(y=Ymatrixsd,x=Xmatrix,standardize=F,intercept=F,family="gaussian",thresh=glmnetthresh)
      lambda_restricted<-matrix(determine_lambdasequence$lambda,nrow=1)
    }
    
    # Time series cross-validation to determine value of the tuning parameter lambda_beta
    cutoff.n<-round(cutoff*nrow(Ymatrix))
    CVscore.beta<-matrix(NA,ncol=ncol(lambda_restricted),nrow=nrow(Ymatrix)-cutoff.n)
    
    for (i in cutoff.n:nrow(Ymatrix)-1){ #Loop to calculate cross-validation score
      # Training Data
      Ytrain<-Ymatrix[1:i,i.r]
      Ytrain.sd<-sd(Ytrain)
      Ytrain.scaled<-Ytrain/Ytrain.sd
      Xtrain<-Xmatrix[1:i,]
      
      # Test Data
      Ytest<-Ymatrix[i+1,i.r]
      Xtest<-Xmatrix[i+1,]
      
      # Estimation
      BETA.scaled<-glmnet(y=Ytrain.scaled,x=Xtrain,lambda=lambda_restricted,standardize=F,intercept=F,family="gaussian",thresh=glmnetthresh)
      B_BETASD<-matrix(BETA.scaled$beta[,which(BETA.scaled$df!=0)],nrow=ncol(Xtrain)) # BETA in standardized scale
      B_BETA<-B_BETASD*Ytrain.sd
      
      B_CV<-apply(B_BETA,2,CVscore.Reg,Y.data=Ytest,X.data=Xtest,Y.sd=Ytrain.sd) # CVscore
      CVscore.beta[i-cutoff.n+1,c(which(BETA.scaled$df!=0))]<-B_CV  
    }
    
    CVscore.AVAILABLE<-as.matrix(CVscore.beta[,apply(CVscore.beta,2,AVAILABLE_LAMBDA)])
    lambda_restricted_AVAILABLE<-lambda_restricted[apply(CVscore.beta,2,AVAILABLE_LAMBDA)]
    lambda.opt<-lambda_restricted_AVAILABLE[which.min(apply(CVscore.AVAILABLE,2,AVERAGE))]
    
    LASSOfinal<-glmnet(y=Ymatrixsd,x=Xmatrix,standardize=F,intercept=F,lambda=lambda.opt,family="gaussian",thresh=glmnetthresh)
    BETAscaled<-matrix(LASSOfinal$beta,ncol=1)
    BETA.Sparse[,i.r]<-BETAscaled*sd(Ymatrix[,i.r])
    
  } # close Loop over cointegration rank
  
  # Determine Omega, conditional on alpha,beta and gamma
  if(intercept==T){
    Resid<-(Y-cbind(1,X)%*%ZBETA)-Z%*%BETA.Sparse%*%t(alpha)
  } else {
    Resid<-(Y-cbind(X)%*%ZBETA)-Z%*%BETA.Sparse%*%t(alpha)
  } 
  
  Resid<-Re(Resid)
  covResid<-cov(Resid)
  if(length(rho.glasso)==1){
    GLASSOfit<-huge(covResid,lambda=rho.glasso,method="glasso",cov.output=T,verbose=F)
    OMEGA<-GLASSOfit$icov[[1]]
  }else{
    GLASSOfit<-huge(Resid,lambda=rho.glasso,method="glasso",cov.output=T,verbose=F)
    huge.BIC<-huge.select(GLASSOfit, criterion = "ebic",verbose=F)
    OMEGA<-huge.BIC$opt.icov
  }

  out<-list(BETA=BETA.Sparse,OMEGA=OMEGA)
}

AVAILABLE_LAMBDA<-function(U){
  # AUXILIARY FUNCTION 
  if (all(is.na(U)==T)) {
    return(F)
  } else {
    return(T)
  }
}

AVERAGE<-function(U){
  # AUXILIARY FUNCTION 
  mean(U[which(is.na(U)==F)])
}

CVscore.Reg<-function(U,X.data,Y.data,Y.sd){
  # AUXILIARY FUNCTION 
  # Standardized Mean Squared  Error
  return(mean((Y.data-X.data%*%U)^2/Y.sd))
}

NORMALIZATION_UNIT<-function(U){
  # AUXILIARY FUNCTION 
  # output: normalized vector U
  length.U<-as.numeric(sqrt(t(U)%*%U))
  if(length.U==0){length.U<-1}
  Uunit<-U/length.U
}

principal_angles<-function(a,b){
  # AUXILIARY FUNCTION 
  # Calculate minimal angle between subspace a and b
  # INPUT: subspace a and b
  # OUTPUT: list of angles
  angles=matrix(0,ncol=ncol(a),nrow=1)
  angles=matrix(0,ncol=ncol(a),nrow=1)
  qa= qr.Q(qr(a))
  qb= qr.Q(qr(b))
  C=svd(t(qa)%*%qb)$d
  rkA=qr(a)$rank;rkB=qr(b)$rank
  if (rkA<=rkB){
    B = qb - qa%*%(t(qa)%*%qb);
  } else {B = qa - qb%*%(t(qb)%*%qa)}
  S=svd(B)$d
  S=sort(S)
  
  for (i in 1:min(rkA,rkB)){
    if (C[i]^2 < 0.5) {angles[1,i]=acos(C[i])}
    else if (S[i]^2 <=0.5) {angles[1,i]=asin(S[i])}
  }
  angles=t(angles)
  
  
  
  return(angles[1])
}

determine_rank<-function(Y,X,Z,r.init=NULL,p,max.iter.lasso=3,conv.lasso=10^-2,max.iter.r=5,beta.init,alpha.init,rho.glasso=0.1,lambda.gamma=matrix(seq(from=0.1,to=0.001,length=10),nrow=1),lambda_beta=matrix(seq(from=2,to=0.001,length=100),nrow=1),glmnetthresh=1e-04){
  ### FUNCTION TO DETERMINE THE COINTEGRATION RANK USING THE RANK SELECTION CRITERION ###
  
  # INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r.init: initial value of cointegration rank
  # p: number of lags to be included
  # max.iter.lasso: maximum number of iterations PML
  # conv.lasso: convergence parameter
  # max.iter.r: maximu number of iterations to compute cointegration rank
  # beta.init: initial value for beta
  # alpha.init: initial value for alpha
  # rho.glasso: tuning parameter for inverse covariance matrix
  # lambda.gamma: tuning parameter for GAMMA
  # lambda.beta: tuning parameter for BETA
  # glmnetthresh: tolerance parameter glmnet function
  
  # OUTPUT
  # rhat: estimated cointegration rank
  # it.r: number of iterations
  # rhat_iterations: estimate of cointegration rank in each iteration
  # mu: value of mu
  # decomp: eigenvalue decomposition
  
  # START CODE
  
  # Starting values
  if(is.null(beta.init)&is.null(alpha.init)){
    q<-ncol(Y)
    beta.init<-matrix(1,ncol=1,nrow=q)
    SVDinit<-svd(t(Z%*%beta.init)%*%(Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T)))
    alpha.init<-t(SVDinit$u%*%t(SVDinit$v))
    beta.init<-matrix(NA,ncol=1,nrow=q)
    Yinit<-(Y-X%*%matrix(rbind(rep(diag(1,ncol(Y)),p-1)),ncol=ncol(Y),byrow=T))%*%alpha.init
    for (i in 1:1){
      FIT<-lars(x=Z,y=Yinit[,i],type="lasso",normalize=F,intercept=F)$beta[-1,]
      beta.init[,i]<-FIT[nrow(FIT),]  #nrow(FIT)
    } 
    beta.init<-apply(beta.init,2,NORMALIZATION_UNIT)
  }
  
  # Preliminaries
  n<-nrow(Y)
  if(is.null(r.init)){
    r.init<-ncol(Y) 
  }
  diff.r<-1
  it.r<-1
  rhat_iterations<-matrix(ncol=max.iter.r,nrow=1)
  
  
  while ((it.r<max.iter.r) & (diff.r>0)){
    # Initialization
    beta.init.fit<-matrix(0,ncol=r.init,nrow=ncol(Y))
    alpha.init.fit<-matrix(0,ncol=r.init,nrow=ncol(Y))
    if (r.init==0){
      beta.init.fit<-matrix(0,ncol=r.init,nrow=ncol(Y))
      alpha.init.fit<-matrix(0,ncol=r.init,nrow=ncol(Y))
      diff.r=0
    } else {
      
      if (ncol(beta.init.fit)>ncol(beta.init)){
        beta.init.fit[,1:ncol(beta.init)]<-beta.init
        alpha.init.fit[,1:ncol(beta.init)]<-alpha.init
      } else {
        beta.init.fit[,1:r.init]<-beta.init[,1:r.init]
        alpha.init.fit[,1:r.init]<-alpha.init[,1:r.init]
      }
      
      # Sparse cointegration fit 
      FITLASSO<-SparseCointegration_RSC(Y=Y,X=X,Z=Z,beta.init=beta.init.fit,alpha.init=alpha.init.fit,p=p,r=r.init,max.iter=max.iter.lasso,conv=conv.lasso,lambda.gamma=lambda.gamma,lambda_beta=lambda_beta,rho.glasso=rho.glasso,glmnetthresh=glmnetthresh)
      
      # Response  and predictors in penalized reduced rank regression
      Y.new<-Y-X%*%FITLASSO$ZBETA
      X.new<-Z
      M<-t(X.new)%*%X.new
      P<-X.new%*%ginv(M)%*%t(X.new)
      l<-ncol(Y.new)
      h<-qr(X.new)$rank
      m<-nrow(Y.new)
      S<-sum((Y.new-P%*%Y.new)^2)/(m*l-h*l)
      mu<-2*S*(l+h)
      decomp<-Re(eigen(t(Y.new)%*%P%*%Y.new)$values)
      khat<-length(which(decomp>mu | decomp==mu ))
      
      # Convergence checking
      rhat_iterations[,it.r]<-khat
      diff.r<-abs(r.init-khat)
      r.init<-khat
      it.r<-it.r+1
    }
  }

  out<-list(rhat=khat,it.r=it.r,rhat_iterations=rhat_iterations,mu=mu,decomp=decomp)
}

SparseCointegration_RSC<-function(p,Y,X,Z,r,alpha.init=NULL,beta.init,max.iter=25,conv=10^-3,lambda.gamma=0.001,lambda_beta=matrix(seq(from=2,to=0.001,length=100),nrow=1),rho.glasso=0.5,cutoff=0.8,intercept=F,glmnetthresh=1e-04){
  ### Sparse cointegration function used in determine_rank function ###
  
  ## INPUT
  # p: number of lagged differences
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r: cointegration rank
  # alpha.init: initial value for adjustment coefficients
  # beta.init: initial value for cointegrating vector
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # lambda.gamma: tuning paramter short-run effects
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # glmnetthresh: tolerance parameter glmnet function
  
  ## OUTPUT
  # BETAhat: estimate of cointegrating vectors 
  # ALPHAhat: estimate of adjustment coefficients
  # ZBETA: estimate of short-run effects
  # OMEGA: estimate of inverse covariance matrix
  
  
  ## START CODE
  
  # Dimensions
  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  #Starting value
  if (is.null(alpha.init)){
    alpha.init<-matrix(rnorm(q*r),ncol=r,nrow=q)
    beta.init<-matrix(rnorm(q*r),ncol=r,nrow=q)
    Pi.init<-alpha.init%*%t(beta.init)
  } else {
    Pi.init<-alpha.init%*%t(beta.init)
  }
  
  # Convergence parameters: initialization
  it<-1
  diff.obj<-10*conv
  Omega.init=diag(1,q)
  value.obj<-matrix(NA,ncol=1,nrow=max.iter+1)
  RESID<-Y-X%*%matrix(rep(diag(1,q),p-1),ncol=q,byrow=T)-Z%*%beta.init%*%t(alpha.init)
  value.obj[1,]<-(1/n)*sum(diag(RESID%*%Omega.init%*%t(RESID)))-log(det(Omega.init))
  
  
  while( (it<max.iter) &  (diff.obj>conv) )
  {
    
    # Obtain Gamma, fixed value of tuning parameter
    FIT1<-NTS.GAMMAFIXED.LassoReg(Y=Y,X=X,Z=Z,Pi=Pi.init,p=p,Omega=Omega.init,lambda.gamma=lambda.gamma,glmnetthresh=glmnetthresh)
    # Obtain alpha
    FIT2<-NTS.ALPHA.Procrusted(Y=Y,X=X,Z=Z,ZBETA=FIT1$ZBETA,r=r,Omega=Omega.init,P=FIT1$P,beta=beta.init)
    # abtain beta and omega
    FIT3<-NTS.BETA(Y=Y,X=X,Z=Z,ZBETA=FIT1$ZBETA,r=r,Omega=Omega.init,P=FIT1$P,alpha=FIT2$ALPHA,alphastar=FIT2$ALPHAstar,lambda_beta=lambda_beta,rho.glasso=rho.glasso,cutoff=cutoff,glmnetthresh=glmnetthresh)
    
    
    # Check convergence
    beta.new<-FIT3$BETA
    if(intercept==T){
      RESID<-Y-X%*%FIT1$ZBETA[-1,]-Z%*%FIT3$BETA%*%t(FIT2$ALPHA)
    }else{
      RESID<-Y-X%*%FIT1$ZBETA-Z%*%FIT3$BETA%*%t(FIT2$ALPHA)
    }
    value.obj[1+it,]<- (1/n)* sum(diag((RESID)%*%FIT3$OMEGA%*%t(RESID))) - log(det(FIT3$OMEGA))
    diff.obj<-abs(value.obj[1+it,]-value.obj[it,])/abs(value.obj[it,])
    beta.init<-matrix(beta.init,nrow=q,ncol=r)
    alpha.init<-FIT2$ALPHA
    beta.init<-beta.new
    Pi.init<-alpha.init%*%t(beta.new) ######Pi.init!!!
    Omega.init<-FIT3$OMEGA
    it<-it+1
  }
  
  out<-list(BETAhat=FIT3$BETA,ALPHAhat=FIT2$ALPHA,it=it,ZBETA=FIT1$ZBETA,OMEGA=FIT3$OMEGA)
}

NTS.GAMMAFIXED.LassoReg<-function(Y,X,Z,Pi,Omega,lambda.gamma=0.001,p,cutoff=0.8,intercept=F,glmnetthresh=1e-04){
  # FUNCTIONS TO ESTIMATE GAMMA, using L1 PENALTY, fixed value of lmabda.gamma
  
  ## INPUT
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # Pi: estimate of cointegration space
  # Omega: estimate of inverse error covariance matrix
  # lambda.gamma: tuning paramter short-run effects
  # p: number of lagged differences
  # cutoff: cutoff value time series cross-validation approach
  # intercept: F do not include intercept, T include intercept
  # glmnetthresh: tolerance parameter glmnet function
  
  ## OUTPUT
  # ZBETA: estimate of short-run effects
  # P: transformation matrix P derived from Omega
  
  
  # Setting dimensions
  q<-ncol(Y)
  
  # Creating data matrices
  Y.matrix<-Y-Z%*%t(Pi)
  X.matrix<-cbind(1,X)
  if(intercept==T){
    X.matrix<-cbind(1,X)
  } else {
    X.matrix<-cbind(X)
  }
  
  
  # Compute transformation matrix P
  decomp<-eigen(Omega)
  P<-(decomp$vectors)%*%diag(sqrt(decomp$values))%*%solve(decomp$vectors)
  
  
  # Final estimate 
  Pfinal<-kronecker(P,diag(rep(1,nrow(X.matrix))))   
  YmatrixP<-Pfinal%*%c(Y.matrix)
  YmatrixP.sd<-sd(Y.matrix)
  YmatrixP.scaled<-YmatrixP/YmatrixP.sd
  XmatrixP<-kronecker(P%*%diag(1,q),diag(rep(1,nrow(X.matrix)))%*%X.matrix)
  
  # Lasso Estimation
  LASSOfinal<-glmnet(y=YmatrixP.scaled,x=XmatrixP,standardize=F,intercept=F,lambda=lambda.gamma,family="gaussian",thresh=glmnetthresh)
  GAMMAscaled<-matrix(LASSOfinal$beta,ncol=1)
  GAMMA.Sparse<-GAMMAscaled*YmatrixP.sd
  
  
  ZBETA<-matrix(NA,ncol=q,nrow=ncol(X.matrix))
  for (i.q in 1:q)
  {
    
    if(intercept==T){
      ZBETA[,i.q]<-GAMMA.Sparse[((i.q-1)*((p-1)*q+1)+1):(i.q*(1+(p-1)*q))]
    } else {
      ZBETA[,i.q]<-GAMMA.Sparse[((i.q-1)*((p-1)*q)+1):(i.q*((p-1)*q))]
    }  
  }  
  
  out<-list(ZBETA=ZBETA,P=P)
}

Bartlettcorrection<-function(r,n,p,k,cancor,S01,V,S11,M02,M22,M12,S00){
  # Reference: Johansen (2002) - Econometrica
  nb<-p-r #number of variables - rank
  
  
  # constants Table I p 1939
  a1_constant<-0.561 # no deterministic trend
  a2_constant<--0.016 # no deterministic trend
  a3_constant<-2.690 # no deterministic trend
  b_constant<--0.569 # no deterministic trend
  
  # 
  # a1_constant<-0.494 # no deterministic trend
  # a2_constant<-0.826 # no deterministic trend
  # a3_constant<-0.829 # no deterministic trend
  # b_constant<--0.200 # no deterministic trend
  
  
  a_expression<- 1 + a1_constant*nb/n + a2_constant*(nb/n)^2 + a3_constant*(nb/n)^3+b_constant/n
  
  ny<-r+(k-1)*p #expression p 1933
  
  
  if (r==0){
    BETAHAT<-matrix(0,ncol=1,nrow=p)
    ALPHAHAT<-matrix(0,ncol=1,nrow=p)
    ALPHAHAT_ORTHOGONAL<-matrix(0,ncol=p-r,nrow=p)
    GAMMA1HAT<-M02%*%solve(M22)
    OMEGAHAT<-S00
    
    P<-GAMMA1HAT
    Q<-diag(1,p)
    
    decomposition<-eigen(P)
    K<-decomposition$vectors
    R<-decomposition$values
    #     round(Re(K%*%diag(R,ncol(P))%*%solve(K)),5)
    
    denumerator<-matrix(NA,ncol=ny,nrow=ny)
    for (i in 1:ny){
      for (j in 1:ny){
        denumerator[i,j]<-Re(1-R[i]*R[j])
      } 
    }
    
    
    
    #     SIGMA<-Re(K%*%Re(solve(K)%*%Q%*%OMEGAHAT%*%t(Q)%*%solve(t(K)))/denumerator%*%t(K))
    SIGMA<-Re(K)%*%(ginv(Re(K))%*%Q%*%OMEGAHAT%*%t(Q)%*%ginv(t(Re(K)))/denumerator)%*%Re(t(K))
    
    
    V_psi<-solve(diag(1,ny)-P)%*%Q%*%OMEGAHAT%*%ALPHAHAT_ORTHOGONAL%*%ginv(t(ALPHAHAT_ORTHOGONAL)%*%OMEGAHAT%*%ALPHAHAT_ORTHOGONAL)%*%t(ALPHAHAT_ORTHOGONAL)%*%OMEGAHAT%*%t(Q)%*%solve(diag(1,ny)-P)%*%ginv(SIGMA)
    c1_constant<-sum(diag(V_psi))
    c2_constant<-sum(diag(diag(1,ny)-solve(diag(1,ny)-P)%*%Q%*%OMEGAHAT%*%t(Q)%*%solve(diag(1,ny)-t(P))%*%ginv(SIGMA)))
    #   c2_constant<-0                                                                                  
    c3_constant<-sum(diag(kronecker((diag(1,ny)-P)%*%V_psi,P)%*%solve(diag(1,ny*ny)-kronecker(P,P))))+ sum(diag(V_psi%*%P%*%solve(diag(1,ny)+P)))
    
    
  } else {
    BETAHAT<-matrix(cancor$ycoef[-(p+1),1:r],ncol=r) #DELETE CONSTANT
    ALPHAHAT<-S01%*%V[,1:r]%*%ginv(t(V[,1:r])%*%S11%*%V[,1:r])
    ALPHAHAT_ORTHOGONAL<- qr.Q(qr(ALPHAHAT),complete=TRUE)[,-r]
    
    GAMMA1HAT<-M02%*%solve(M22)-ALPHAHAT%*%t(V[,(1:r)])%*%M12%*%solve(M22)
    OMEGAHAT<-S00-ALPHAHAT%*%(t(BETAHAT)%*%S11[-(p+1),-(p+1)]%*%BETAHAT)%*%t(ALPHAHAT) #DELETE CONSTANT
    
    P<- rbind(cbind((diag(1,r)+t(BETAHAT)%*%ALPHAHAT),t(BETAHAT)%*%GAMMA1HAT),cbind(ALPHAHAT,GAMMA1HAT))
    Q<-rbind(t(BETAHAT),diag(1,p))
    
    decomposition<-eigen(P)
    K<-decomposition$vectors
    R<-decomposition$values
    #     round(Re(K%*%diag(R,ncol(P))%*%solve(K)),5)
    
    
    denumerator<-matrix(NA,ncol=ny,nrow=ny)
    for (i in 1:ny){
      for (j in 1:ny){
        denumerator[i,j]<-Re(1-R[i]*R[j])
      }
      
    }
    
    #     SIGMA<-Re(K%*%Re(solve(K)%*%Q%*%OMEGAHAT%*%t(Q)%*%solve(t(K)))/denumerator%*%t(K))
    SIGMA<-Re(K)%*%(ginv(Re(K))%*%Q%*%OMEGAHAT%*%t(Q)%*%ginv(t(Re(K)))/denumerator)%*%Re(t(K))
    V_psi<-solve(diag(1,ny)-P)%*%Q%*%OMEGAHAT%*%ALPHAHAT_ORTHOGONAL%*%ginv(t(ALPHAHAT_ORTHOGONAL)%*%OMEGAHAT%*%ALPHAHAT_ORTHOGONAL)%*%t(ALPHAHAT_ORTHOGONAL)%*%OMEGAHAT%*%t(Q)%*%solve(diag(1,ny)-P)%*%ginv(SIGMA)
    c1_constant<-sum(diag(V_psi))
    c2_constant<-sum(diag(diag(1,ny)-solve(diag(1,ny)-P)%*%Q%*%OMEGAHAT%*%t(Q)%*%solve(diag(1,ny)-t(P))%*%ginv(SIGMA)))
    c3_constant<-sum(diag(kronecker((diag(1,ny)-P)%*%V_psi,P)%*%solve(diag(1,ny*ny)-kronecker(P,P))))+ sum(diag(V_psi%*%P%*%solve(diag(1,ny)+P)))
    
  }
  h1_constant<-0
  h2_constant<-0
  h3_constant<-0
  g0_constant<--0.506
  g1_constant<-0.02
  g2_constant<-0.07
  g3_constant<--0.144
  
  # h1_constant<-0
  # h2_constant<-0.197
  # h3_constant<-0.036
  # g0_constant<--0.496
  # g1_constant<-0.166
  # g2_constant<-0.079
  # g3_constant<--0.076
  
  
  h_expression<-h1_constant/nb + h2_constant/(nb^2)+h3_constant/(nb^3)
  g_expression<-g0_constant+g1_constant/nb + g2_constant/(nb^2)+g3_constant/(nb^3)
  # b_expression<-c1_constant*(1+h_expression)+ (nb*c2_constant+2*(c3_constant))*g_expression/(nb^2) #check nd
  b_expression<-c1_constant*(1+h_expression)+ (nb*c2_constant+2*(c3_constant+nb*c1_constant))*g_expression/(nb^2) #check nd
  
  correctionterm<-a_expression*(1+b_expression/n)
  
  out<-list(a_expression=a_expression,b_expression=b_expression,correctionterm=correctionterm)
}

BOOTSTRAP<-function(r,V,p,M02,M22,M12,LEVEL_X_AUX,DESIGN_X_AUX,DELTA_X_AUX,Nboot,test_stat_Johansen){
  ### Function: Bootstrap procedure of Cavaliere et al. (2012) to determine cointegration rank ###
  
  ### START CODE
  NOTOK<-0
  if (r==0){
    BETA_r<-matrix(0,ncol=1,nrow=p)
    ALPHA_r<-matrix(0,ncol=1,nrow=p)
    RHO_r<-0
    GAMMA1_HAT<-M02%*%solve(M22)
    
    BETA_FIT<-rbind(t(ALPHA_r%*%t(BETA_r)),t(GAMMA1_HAT),t(ALPHA_r*RHO_r))
    BIG_X<-cbind(LEVEL_X_AUX[,(1:p)],DESIGN_X_AUX,LEVEL_X_AUX[,ncol(LEVEL_X_AUX)])
    DELTA_FIT<-BIG_X%*%BETA_FIT
    resid_HAT<-DELTA_X_AUX-DELTA_FIT
  } else { 
    BETA_r<-V[1:p,(1:r)]
    RHO_r<-matrix(V[nrow(V),(1:r)],ncol=r)
    ALPHA_r<-S01%*%V[,(1:r)]%*%solve(t(V[,(1:r)])%*%S11%*%V[,(1:r)])
    
    GAMMA1_HAT<-M02%*%solve(M22)-ALPHA_r%*%t(V[,(1:r)])%*%M12%*%solve(M22)
    
    BETA_FIT<-rbind(t(ALPHA_r%*%t(BETA_r)),t(GAMMA1_HAT),t(ALPHA_r%*%t(RHO_r)))
    BIG_X<-cbind(LEVEL_X_AUX[,(1:p)],DESIGN_X_AUX,LEVEL_X_AUX[,ncol(LEVEL_X_AUX)])
    DELTA_FIT<-BIG_X%*%BETA_FIT
    resid_HAT<-DELTA_X_AUX-DELTA_FIT  
  }
  BOOT_test_stat_r<-matrix(data=0,ncol=2,nrow=Nboot)
  
  A1j<-diag(1,p)+ALPHA_r%*%t(BETA_r)+GAMMA1_HAT
  A2<--GAMMA1_HAT
  companion<-cbind(rbind(A1j,diag(1,p)),rbind(A2,matrix(0,nrow=p,ncol=p)))
  eigenvalues<-round(abs(eigen(companion)$values),digits=1)
  eigenvalues_OK<-eigenvalues<=1
  if (sum(eigenvalues_OK)!=length(eigenvalues_OK)){NOTOK=NOTOK+1}
  
  for (nboot in 1:Nboot){
    BootIndex <- ceiling(runif(n,min=0,max=nrow(resid_HAT)))
    BootResiduals <- resid_HAT[BootIndex,]
    
    BOOT_X<-matrix(0,ncol=p,nrow=n)
    BOOT_X[1:2,]<-0
    for (i in 3:n){
      BOOT_X[i,]<- matrix(data=BOOT_X[i-1,],ncol=1) + ALPHA_r%*%t(BETA_r)%*%matrix(data=BOOT_X[i-1,],ncol=1) + GAMMA1_HAT%*%matrix(data=(BOOT_X[i-1,]-BOOT_X[i-2,]),ncol=1) + BootResiduals[i,]
    }
    
    BOOT_X<-stdize(BOOT_X)
    
    BOOT_DELTA_X<-diff(BOOT_X,differences=1)
    BOOT_LEVEL_X<-BOOT_X[-nrow(BOOT_X),]
    BOOT_data<-embed(BOOT_DELTA_X,dimension=k)
    
    BOOT_DELTA_X_AUX<-BOOT_data[,(1:p)]
    BOOT_LEVEL_X_AUX<-cbind(BOOT_LEVEL_X[-(1:(k-1)),],1)
    BOOT_DESIGN_X_AUX<-cbind(BOOT_data[,-(1:p)])
    
    BOOT_BETA_DELTA<-solve(t(BOOT_DESIGN_X_AUX)%*%BOOT_DESIGN_X_AUX)%*%t(BOOT_DESIGN_X_AUX)%*%BOOT_DELTA_X_AUX
    BOOT_BETA_LEVEL<-solve(t(BOOT_DESIGN_X_AUX)%*%BOOT_DESIGN_X_AUX)%*%t(BOOT_DESIGN_X_AUX)%*%BOOT_LEVEL_X_AUX
    
    BOOT_RESID_DELTA<-BOOT_DELTA_X_AUX-BOOT_DESIGN_X_AUX%*%BOOT_BETA_DELTA
    BOOT_RESID_LEVEL<-BOOT_LEVEL_X_AUX-BOOT_DESIGN_X_AUX%*%BOOT_BETA_LEVEL
    
    # Canonical Correlation Analysis
    BOOT_cancor=cancor(BOOT_RESID_DELTA,BOOT_RESID_LEVEL)
    BOOT_eigenvalues=BOOT_cancor$cor^2
    
    BOOT_test_stat_r[nboot,1]<-sum(-n*log(1-BOOT_eigenvalues[(r+1):p]))
    
    if (BOOT_test_stat_r[nboot,1] > test_stat_Johansen[r+1,1] ) {BOOT_test_stat_r[nboot,2]=1}
  }
  BOOT_fraction_r<-sum(BOOT_test_stat_r[,2])/Nboot
  out=list(detail=BOOT_test_stat_r,fraction=BOOT_fraction_r,Bootstrap_flag=NOTOK)
}
