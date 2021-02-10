#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' @rdname systematic_resampling
//' @title Perform systematic resampling
//' @param normweights normalized particle weights 
//' @param nparticles number of particles
//' @param uniform uniform variable between 0 to 1 / nparticles 
//' @return \code{ancestors} ancestor indices
//' @export
// [[Rcpp::export]]
arma::vec systematic_resampling(const arma::vec & normweights, int nparticles, double uniform) {
  arma::vec ancestors = arma::zeros(nparticles);
  int j = 0;
  double csw = normweights(0);
  for(int k = 0; k < nparticles; k++){
    while(csw < uniform){
      j++;
      csw += normweights(j);
    }
    uniform = uniform + 1. / nparticles;
    ancestors(k) = j + 1;
  }
  return ancestors;
}

