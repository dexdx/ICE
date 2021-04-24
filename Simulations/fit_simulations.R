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
#' days to execute entirely, in particular to the sPSC. Summarized results of these
#' fitted processes are available on https://jonkq.github.io/research.html, under
#' Master Dissertation.
#' 
##########################################################################
##### FIT + CROSS VALIDATION #############################################
##########################################################################
source('sVAR estimators.R')
X.list = list(X1.list,X2.list,X3.list,X4.list)
B.list = list(B1.full,B2.full,B3.full,B4.full)