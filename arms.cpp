# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

extern "C" {
#include "arms.c"
}

using namespace Rcpp;
using namespace arma;
using namespace std;

struct log_p_param {
  vec theta_i;
  mat covmat_inv;
  unsigned int k;
};

double log_den_func(double x, void *params){
  log_p_param *d;
  d = static_cast<log_p_param *> (params);

  mat covmat_inv = d->covmat_inv;
  unsigned int k = d->k;
  vec theta_i = d->theta_i;
  theta_i(k) = x;
  
  double logden = - 0.5 * as_scalar(theta_i.t() * covmat_inv * theta_i);
  return(logden);
}

/* *********************************************************************** */

// [[Rcpp::export]]
mat sampling(mat covmat, int N){
  
  int K = covmat.n_cols;
  mat theta = zeros(N, K);
  mat covmat_inv = inv_sympd(covmat);
  
  int err, ninit = 4, npoint = 100, nsamp = 1, ncent = 0;
  int neval;
  double xinit[10]={-3.0, -1.0, 1.0, 3.0}, xl = -50.0, xr = 50.0;
  double xsamp[100], xcent[10], qcent[10] = {5., 30., 70., 95.};
  double convex = 1.;
  int dometrop = 0;
  // Set dometrop=0 if you know the log-density is concave; 
  // if you know it isn't, or you're not sure, set dometrop=1.
  double xprev = 0.0;
  
  /* initialise random number generator */
  // unsigned seed = 44;
  // srand(seed);
  
  log_p_param params_data;
  params_data.covmat_inv = covmat_inv;
  params_data.theta_i = zeros(K);
  
  for(int i=0; i<N; ++i){
    for(int k=0; k<K; ++k){
      params_data.k = k;
      err = arms(xinit,ninit,&xl,&xr,log_den_func,&params_data,&convex,
                 npoint,dometrop,&xprev,xsamp,nsamp,qcent,xcent,ncent,&neval);
      if(err>0){
        Rprintf("error code: %d", err);
        Rcpp::stop("\n");
      }
      params_data.theta_i(k) = xsamp[0];
    }
    theta.row(i) = params_data.theta_i.t();
  }
  
  return (theta);
}









