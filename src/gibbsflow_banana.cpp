#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Problem settings
const int dimension = 2;
const double prior_sd = 1.0;

// Prior 
const double log2pi = std::log(2.0 * arma::datum::pi);
const arma::rowvec prior_mean = arma::zeros<arma::rowvec>(dimension);
const arma::mat prior_cov = prior_sd * arma::eye(dimension, dimension);
const arma::mat prior_precision = arma::inv_sympd(prior_cov);
const arma::mat prior_chol = arma::chol(prior_cov);
const arma::mat prior_rooti = arma::trans(arma::inv(arma::trimatu(prior_chol)));
const double prior_rootisum = arma::sum(log(prior_rooti.diag()));

//' @rdname banana_logprior
//' @title Evaluate banana example prior density
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec banana_logprior(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = -(static_cast<double>(dimension)/2.0) * log2pi + prior_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = prior_rooti * arma::trans(x.row(i) - prior_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname banana_gradlogprior
//' @title Evaluate gradient of banana example prior density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat banana_gradlogprior(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (prior_mean - x.row(i)) * prior_precision;
  }
  
  return(output);
}

//' @rdname banana_partialderivative_logprior
//' @title Evaluate partial derivative of banana example prior density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double banana_partialderivative_logprior(int dim, arma::rowvec x){
  arma::rowvec xc = prior_mean - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += prior_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname banana_sampleprior
//' @title Sample from banana example prior distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat banana_sampleprior(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = prior_mean + Z * prior_chol;
  }
  
  return(output);
}

// Likelihood
//' @rdname banana_loglikelihood
//' @title Evaluate banana example loglikelihood function
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec banana_loglikelihood(arma::mat x, double like_alpha, double like_beta){
  int n = x.n_rows;
  arma::vec output(n);
  double term1, term2;
  
  for (int i = 0; i < n; i++){
    term1 = like_alpha - x(i, 0);
    term2 = x(i, 1) - x(i, 0) * x(i, 0);
    output(i) = - term1 * term1 - like_beta * term2 * term2;
  }
  
  return(output);
}

//' @rdname banana_gradloglikelihood
//' @title Evaluate gradient of banana example loglikelihood function
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat banana_gradloglikelihood(arma::mat x, double like_alpha, double like_beta){
  int n = x.n_rows;
  double term1, term2;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    term1 = like_alpha - x(i, 0);
    term2 = x(i, 1) - x(i, 0) * x(i, 0);
    output(i, 0) = 2.0 * term1 + 4.0 * like_beta * x(i, 0) * term2;
    output(i, 1) = - 2.0 * like_beta * term2;
  }
  
  return(output);
}

//' @rdname banana_partialderivative_loglikelihood
//' @title Evaluate partial derivative of banana example likelihood function
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double banana_partialderivative_loglikelihood(int dim, arma::rowvec x, double like_alpha, double like_beta){
  double output;
  
  switch (dim) {
  case 0 : 
    output = 2.0 * ( like_alpha - x(0) ) + 4.0 * like_beta * x(0) * ( x(1) - x(0) * x(0) );
    break;
  case 1 : 
    output = - 2.0 * like_beta * ( x(1) - x(0) * x(0) );
    break;
  }
  
  return(output);
  
}

// Gibbs flow computation
const double lowerbound = - 6.0;
const double upperbound = 10.0;
const int ngridpoints = 200;
const arma::vec grid = arma::linspace(lowerbound, upperbound, ngridpoints);
const double mesh = grid(1) - grid(0);

//' @rdname banana_gibbsflow
//' @title Compute Gibbs flow for banana example
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List banana_gibbsflow(double stepsize, double lambda, double derivative_lambda,
                                    arma::mat xparticles, arma::vec logdensity, double like_alpha, double like_beta){
  int N = xparticles.n_rows; // no. of particles

  // Output
  arma::rowvec log_jacobian_dets(N);

  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double log_jacobian_det, loglikelihood_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative;
  arma::rowvec xpoint, gridpoint;
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;

  for (int n = 0; n < N; n++){

    xpoint = xparticles.row(n);
    log_jacobian_det = 0;

    for (int i = 0; i < dimension; i++){
      loglikelihood_x = as_scalar(banana_loglikelihood(xpoint, like_alpha, like_beta));
      target_x = std::exp( as_scalar(banana_logprior(xpoint)) + lambda * loglikelihood_x );
      loglikelihood_target_x = loglikelihood_x * target_x;

      // Construct quadrature grids
      xindex = static_cast<int>((xpoint(i) - lowerbound) / mesh) + 1;
      nlowergrid = xindex + 1;
      lowergrid = arma::zeros<arma::vec>(nlowergrid);
      lowergrid(arma::span(0, xindex - 1)) = grid(arma::span(0, xindex - 1));
      lowergrid(xindex) = xpoint(i);

      nuppergrid = ngridpoints - xindex + 1;
      uppergrid = arma::zeros<arma::vec>(nuppergrid);
      uppergrid(0) = xpoint(i);
      uppergrid(arma::span(1, nuppergrid - 1)) = grid(arma::span(xindex, ngridpoints - 1));

      gridpoint = xpoint;

      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);

      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;

      for (int p = 0; p < xindex; p++){
        gridpoint(i) = lowergrid(p);
        loglikelihood_gridpoint = as_scalar(banana_loglikelihood(gridpoint, like_alpha, like_beta));
        target_gridpoint = std::exp( as_scalar(banana_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = loglikelihood_gridpoint * target_gridpoint;
      }

      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);

      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint(i) = uppergrid(p);
        loglikelihood_gridpoint = as_scalar(banana_loglikelihood(gridpoint, like_alpha, like_beta));
        target_gridpoint = std::exp( as_scalar(banana_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
        target_uppergrid(p) = target_gridpoint;
        loglikelihood_target_uppergrid(p) = loglikelihood_gridpoint * target_gridpoint;
      }

      // Compute conditional CDF
      lowerintegral_target = as_scalar(trapz(lowergrid, target_lowergrid));
      upperintegral_target = as_scalar(trapz(uppergrid, target_uppergrid));
      integral_target = lowerintegral_target + upperintegral_target;
      conditional_cdf = lowerintegral_target / integral_target;

      // Compute expected loglikelihood
      lowerintegral_loglikelihood = as_scalar(trapz(lowergrid, loglikelihood_target_lowergrid));
      upperintegral_loglikelihood = as_scalar(trapz(uppergrid, loglikelihood_target_uppergrid));
      integral_loglikelihood = lowerintegral_loglikelihood + upperintegral_loglikelihood;

      // Compute Gibbs velocity field and its partial derivative
      gibbs_velocity = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
      logtarget_partial_derivative = banana_partialderivative_logprior(i, xpoint) + lambda * banana_partialderivative_loglikelihood(i, xpoint, like_alpha, like_beta);
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - loglikelihood_x) -
        gibbs_velocity * logtarget_partial_derivative;
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
      xpoint(i) = xpoint(i) + stepsize * gibbs_velocity;

    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;

  }

  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles,
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);

}

