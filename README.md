# ARMScore
## Stand alone version of ARMS sampler

Generated from the C code of original ARMS algorithm by [Wally Gilks](https://www1.maths.leeds.ac.uk/~wally.gilks/adaptive.rejection/web_page/Welcome.html).

## Example for sampling multivariate normal
```
library(Rcpp)
library(RcppArmadillo)

sourceCpp('arms.cpp')

covmat <- matrix(0.2, 4, 4)
diag(covmat) <- c(1,2,3,4)
covmat

out = sampling(covmat, 1e5)
cov(out[-(1:2000),])

library(mvtnorm)
out2 = rmvnorm(1e5, sigma = covmat)
cov(out2)
```

## How to use
- Define sampling procedure (possibly with C++ `for` loop) in `arms.cpp` file
- Invoke exposed `sampling()` function in R
- Modify `arms.cpp` to sample distributions through MCMC (Gibbs).
- Check `arms.txt` as a document about `arms()` function arguments given by the original author.
