#' The following estimators will be fitted to all of the simulated VAR(1) paths
#' with the exception of the Johansen's unregularized estimator (4) that will
#' only be fitted to the low-dimensional set-ups since it is unfeasible in 
#' high dimensions.
#' 
#' ICE:
#' 1 : MULTIVARIATE LASSO
#' 2 : SPARSE YULE-WALKER (sYWE)
#' 3 : SPARSE SPECTRAL COHERENCE (sPSC)
#' DiCE:
#' 4 : JOHANSEN (DiCE)
#' 5 : SPARSE JOHANSEN (sDiCE)
#' 
#' For each fit, the L1 norm of the difference between the estimated cointegrating
#' vector and the true cointegrating vector will be saved for performance comparison
#' between estimators across the different set-ups.
#' 
#' In order not to save the 400 simulated processes, the 'processes_simulations.R'
#' file must be run before running this file, so that the processes lie in the 
#' current environment and are ready to be used.
#' 
#' Warning: This file's primary purpose is educational. Since some of the regularized
#' estimators are computationally demanding, running this file will take several
#' days to execute entirely, in particular due to the sPSC. Summarized results of
#' these fitted processes are available on https://jonkq.github.io/research.html,
#' under 'Master Dissertation'.


##########################################################################
##### FIT + CROSS VALIDATION #############################################
##########################################################################
source('sVAR estimators.R')
source('ICE.R')
X.list = list(X1.list,X2.list,X3.list,X4.list)
B.list = list(B1.full,B2.full,B3.full,B4.full)


##### 1) LASSO ###########################################################

# initliaze performance matrix of LASSO
L1B.lasso = matrix(0,M,length(X.list))

for(l in 1:length(X.list)){
  for(m in 1:M){
    
    x = X.list[[l]][[m]]
    N = dim(x)[1]
    K = dim(x)[2]
    x.train = x[1:(N/2),]    # fitting set
    x.test  = x[(N/2+1):N,]  # (cross) validation set
    
    # demean sets with estimate of train
    train.m = colMeans(x.train)
    x.train.dm = t(apply(x.train, 1, function(x) x - train.m))
    test.m = colMeans(x.test)
    x.test.dm = t(apply(x.test, 1, function(x) x - test.m))
    
    A.tune = MVlasso(x.train.dm) # estimated transition matrices as a fct of lasso tuning parameter
    num.tune = length(A.tune)
    x.fitted.tune = vector(mode = 'list', length = num.tune)
    resid.tune = vector(mode = 'list', length = num.tune)
    ssr.tune = numeric(num.tune)
    
    
    ## Two 2-fold cross-validation (cv) methods to choose optimal A.tune
    ## Comment out either one
    
    #1: 2-fold cv BY FORECAST
    for(u in 1:num.tune){
      x.fitted.tune[[u]] = VAR1.fit(x.test[1,], nrow(x.test), A.tune[[u]])
      resid.tune[[u]] = x.test - x.fitted.tune[[u]]
    }
    ssr.tune = sapply(resid.tune, function(x) sum(x^2))
    ssr.ord = order(ssr.tune)

    if(ssr.ord[1] == 1){
      s = 2
      while(ssr.ord[s] == ssr.ord[s-1] + 1 && s < num.tune) s = s + 1
      tune.min = ssr.ord[s]
    } else tune.min = which.min(ssr.tune)
    
    
    # #2: 2-fold cv BY SIMULATION
    # S = 50
    # ssr.sim = matrix(0,nrow=S,ncol=num.tune)
    # #ssr.avg = numeric(num.tune)
    # ssr.med = numeric(num.tune)
    # 
    # for(u in 1:num.tune){
    #   for(s in 1:S){
    #     x.fitted.tune = VARMAsim(nrow(x.test), arlags = 1, phi = A.tune[[u]], skip = 0, sigma = diag(1,K))$series
    #     x.fitted.tune = t(apply(x.fitted.tune, 1, function(x) x + x.test[1,]))
    #     resid.sim = x.test - x.fitted.tune
    #     ssr.sim[s,u] = sum(resid.sim^2)
    #   }
    #   #ssr.avg[u] = mean(ssr.sim[,u])
    #   ssr.med[u] = median(ssr.sim[,u])
    # }
    # 
    # if(l %in% c(3,4)){ #Sparse structure
    #   tune.min = which.min(ssr.med[1:40])  
    # } else{ #Non-Sparse structure
    #   tune.min = which.min(ssr.med[41:num.tune]) + 41 - 1
    # }
    
    
    ## Retrieve estimated cointegrating vector with ICE
    A = A.tune[[tune.min]]
    PI = A - diag(1,K)
    B = ICE(PI, 1)

    
    ## Save L1 norm of the difference between true and estimated cointegrating vector
    L1B.lasso[m,l] = sum(abs(B.list[[l]][,m] - B))
    write.table(L1B.lasso, file = 'Output/performance results/L1B.lasso.txt')

  }
}



##### 2) sYWE ############################################################
L1B.sYWE = matrix(0,M,length(X.list))


for(l in 1:length(X.list)){
  for(m in 1:M){
    
    start = proc.time()[3]
    end = start + 60
    done = FALSE
    
    while(proc.time()[3] < end && done == FALSE){
      
      x = X.list[[l]][[m]]
      N = dim(x)[1]
      K = dim(x)[2]
      x.train = x[1:(N/2),]
      x.test  = x[(N/2+1):N,]
      
      train.m = colMeans(x.train)
      x.train.dm = t(apply(x.train, 1, function(x) x - train.m))
      
      A.tune = sYWE(x.train.dm)
      
      num.tune = length(A.tune)
      x.fitted.tune = vector(mode = 'list', length = num.tune)
      resid.tune = vector(mode = 'list', length = num.tune)
      ssr.tune = numeric(num.tune)
      
      for(u in 1:num.tune){
        x.fitted.tune[[u]] = VAR1.fit(x.test[1,], nrow(x.test), A.tune[[u]])
        resid.tune[[u]] = x.test - x.fitted.tune[[u]]
      }
      
      ssr.tune = sapply(resid.tune, function(x) sum(x^2))
      ssr.ord = order(ssr.tune, decreasing = T)
      if(ssr.ord[1] == 1){
        s = 2
        while(ssr.ord[s] == ssr.ord[s-1] + 1 && s < num.tune) s = s + 1
        tune.min = ssr.ord[s]
      } else tune.min = which.min(ssr.tune)
      

      A = A.tune[[tune.min]]
      PI = A - diag(1,K)
      B = ICE(PI, 1)
      
      L1B.sYWE[m,l] = sum(abs(B.list[[l]][,m] - B))
      write.table(L1B.sYWE, file = 'Output/performance results/L1B.sYWE.txt')
      
      done = TRUE
      
    }
  }
}


##### 3) sPSC ############################################################
L1B.sPSC = matrix(0,M,length(X.list))

for(l in 1:length(X.list)){
  for(m in 1:M){
    
    start = proc.time()[3]
    end = start + 600
    done = FALSE
    
    while(proc.time()[3] < end && done == FALSE){
      
      x = X.list[[l]][[m]]
      N = dim(x)[1]
      K = dim(x)[2]
      x.train = x[1:(N/2),]
      x.test  = x[(N/2+1):N,]
      
      A.tune = sPSC(x.train)
      
      num.tune = length(A.tune)
      x.fitted.tune = vector(mode = 'list', length = num.tune)
      resid.tune = vector(mode = 'list', length = num.tune)
      ssr.tune = numeric(num.tune)
      
      for(u in 1:num.tune){
        
        x.fitted.tune[[u]] = VAR1.fit(x.test[1,], nrow(x.test), A.tune[[u]])
        resid.tune[[u]] = x.test - x.fitted.tune[[u]]
        ssr.tune = sapply(resid.tune, function(x) sum(x^2))
        ssr.ord = order(ssr.tune)
        
        if(ssr.ord[1] == 1){
          s = 2
          while(ssr.ord[s] == ssr.ord[s-1] + 1 && s < num.tune) s = s + 1
          tune.min = ssr.ord[s]
        } else tune.min = which.min(ssr.tune)
        
      }
      
      A = A.tune[[tune.min]]
      PI = A - diag(1,K)
      B = ICE(PI, 1)

      L1B.sPSC[m,l] = sum((B.list[[l]][,m] - B)^2)
      write.table(L1B.sPSC, file = 'Output/performance results/L1B.sPSC.txt')
      
      done = TRUE
      
    }
  }
}


##### 4) DiCE ############################################################
L1B.dice = matrix(0,M,length(X.list))

for(l in 1:length(X.list)){
  for(m in 1:M){
    
    x = X.list[[l]][[m]]
    N = dim(x)[1]
    K = dim(x)[2]
    # x.train = x[1:(N/2),]
    # x.test  = x[(N/2+1):N,]
    
    dice = ca.jo(x, type='trace', ecdet='none', spec='transitory')
    B = dice@V[,1]
    L1B.dice[m,l] = sum((B.list[[l]][,m] - B)^2)
    write.table(L1B.dice, file = 'Output/performance results/L1B.dice.txt')
    
  }
}


##### 5) sDiCE ###########################################################
L1B.sdice = matrix(0,M,length(X.list))
source('Wilms_Functions.R')

for(l in 1:length(X.list)){
  for(m in 1:M){
    
    x = X.list[[l]][[m]]
    x = x[1:floor(nrow(x)/2),]
    
    N = dim(x)[1]
    K = dim(x)[2]
    p = 2
    r = 1
    
    dsY = embed(diff(x),p)
    Z = x[p:(N-1),]
    Y = dsY[,1:K]
    X = dsY[,-(1:K)]
    
    dice = try(SparseCointegration_Lasso(Y=Y,X=X,Z=Z,beta.init=NULL,alpha.init=NULL,p=p,r=r,max.iter=10,conv=10^-2,lambda.gamma=matrix(seq(from=0.01,to=0.001,length=5),nrow=1),rho.glasso=seq(from=1,to=0.1,length=5),lambda_beta=matrix(seq(from=0.1,to=0.001,length=100),nrow=1)))
    if(class(dice) == 'try-error'){
      L1B.dice.mat[m,l] = NA
      erate.dice.mat[m,l] = NA
      write.table(L1B.dice.mat, file = 'Output/performance results/L1B.sDICE.txt')
    } else{
      B = dice$BETAhat
      L1B.dice.mat[m,l] = sum((B.list[[l]][,m] - B)^2)
      write.table(L1B.dice.mat, file = 'Output/performance results/L1B.sDICE.txt')
    }
    
  }
}
