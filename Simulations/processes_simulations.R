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
