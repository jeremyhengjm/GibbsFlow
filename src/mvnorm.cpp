#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

const double log2pi = std::log(2.0 * arma::datum::pi);

//' @rdname mvnpdf
//' @title Evaluate multivariate Gaussian density
//' @param x evaluation points  
//' @param mu mean vector  
//' @param sigma covariance matrix
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mvnpdf(arma::mat x, arma::rowvec mu, arma::mat sigma){
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::vec output(n);
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(dim)/2.0) * log2pi + rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = rooti * arma::trans(x.row(i) - mu);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname mvn_cholesky_factor
//' @title Compute Cholesky factor
//' @param mu mean vector  
//' @param sigma covariance matrix
//' @return density values
//' @export
// [[Rcpp::export]]
Rcpp::List mvn_cholesky_factor(arma::mat sigma){
  arma::mat cholesky = arma::chol(sigma);
  arma::mat rooti = arma::trans(arma::inv(arma::trimatu(cholesky)));
  double rootisum = arma::sum(log(rooti.diag()));
  
  return Rcpp::List::create(Rcpp::Named("chol") = cholesky,
                            Rcpp::Named("rooti") = rooti,
                            Rcpp::Named("rootisum") = rootisum);  
  
}

//' @rdname mvnpdf_chol
//' @title Evaluate multivariate Gaussian density (with pre-computed Cholesky factor)
//' @param x evaluation points  
//' @param mu mean vector  
//' @param rooti inverse of Cholesky factor
//' @param rootisum normalizing constant term
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mvnpdf_chol(arma::mat x, arma::rowvec mu, arma::mat rooti, double rootisum){
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::vec output(n);
  double constants = -(static_cast<double>(dim)/2.0) * log2pi + rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = rooti * arma::trans(x.row(i) - mu);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname gradmvnpdf
//' @title Evaluate gradient of multivariate Gaussian prior density
//' @param x evaluation points  
//' @param mu mean vector  
//' @param precision precision matrix
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat gradmvnpdf(arma::mat x, arma::rowvec mu, arma::mat precision){
  int n = x.n_rows;
  int dim = x.n_cols;
  arma::mat output(n, dim);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (mu - x.row(i)) * precision;
  }
  
  return(output);
}

//' @rdname mvnrnd
//' @title Simulate from multivariate Gaussian distribution
//' @param n number of samples
//' @param mu mean vector  
//' @param sigma covariance matrix
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat mvnrnd(int n, arma::rowvec mu, arma::mat sigma){
  int dim = mu.n_elem;
  arma::mat Z = arma::randn(n, dim);
  return(arma::repmat(mu,n,1) + Z * arma::chol(sigma));
}

//' @rdname mvnrnd_chol
//' @title Simulate from multivariate Gaussian distribution (with pre-computed Cholesky factor)
//' @param n number of samples
//' @param mu mean vector  
//' @param chol Cholesky factor of covariance matrix
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat mvnrnd_chol(int n, arma::rowvec mu, arma::mat chol){
  int dim = mu.n_elem;
  arma::mat Z = arma::randn(n, dim);
  return(arma::repmat(mu,n,1) + Z * chol);
}
