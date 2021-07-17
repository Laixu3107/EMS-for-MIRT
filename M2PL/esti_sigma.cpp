#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double obj_func_cpp(arma::mat sigma, arma::mat sigma_hat){
  arma::mat sigma_inv = arma::inv(sigma);
  return arma::accu( sigma_inv % sigma_hat ) + log(arma::det(sigma));
}

// [[Rcpp::export]]
arma::mat calcu_sigma_cmle_cpp(arma::mat scov, arma::mat scor, double tol){

  arma::mat sigma_hat = scov;
  arma::mat sigma0 = scor;
  arma::mat sigma1 = sigma0;
  arma::mat tmp = sigma0;
  double eps = 1;
  double step = 1;
  while(eps > tol){
    step = 1;
    tmp = arma::inv(sigma0);
    sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
    sigma1.diag().ones();
    sigma1 = arma::symmatu(sigma1);   // add 2021.04.25
    while(obj_func_cpp(sigma0, sigma_hat) < obj_func_cpp(sigma1, sigma_hat) ||
          min(arma::eig_sym(sigma1)) < 0){
      step *= 0.5;
      sigma1 = sigma0 - step * ( - tmp * sigma_hat * tmp + tmp );
      sigma1.diag().ones();
      sigma1 = arma::symmatu(sigma1);   // add 2021.04.25
    }
    eps = obj_func_cpp(sigma0, sigma_hat) - obj_func_cpp(sigma1, sigma_hat);
    // Rprintf("eps= %f\n", eps);
    sigma0 = sigma1;
  }
  return sigma0;
}

