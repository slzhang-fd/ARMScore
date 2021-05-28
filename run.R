setwd("~/Downloads/testarms")

library(Rcpp)
library(RcppArmadillo)

sourceCpp('arms.cpp')

covmat <- matrix(0.4, 4, 4)
diag(covmat) <- c(1,2,3,4)
covmat

out = sampling(covmat, 1e5)
cov(out[-(1:2000),])

library(mvtnorm)
out2 = rmvnorm(1e5, sigma = covmat)
cov(out2)


