#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Problem settings
const int dimension = 2;
const int ncomponents = 4;
const double like_xi = 6.0;
const double like_rho = 0.75;
const double like_var = 1.0;
arma::rowvec construct_like_weights(int ncomponents){
  arma::rowvec output(ncomponents);
  output(0) = 0.4;
  output(1) = 0.1;
  output(2) = 0.4;
  output(3) = 0.1;
  return (output);
  
}
const arma::rowvec like_weights = construct_like_weights(ncomponents);

// Prior 
const double log2pi = std::log(2.0 * arma::datum::pi);
const arma::rowvec prior_mean = arma::zeros<arma::rowvec>(dimension);
const arma::mat prior_cov = arma::eye(dimension, dimension);
const arma::mat prior_precision = arma::inv_sympd(prior_cov);
const arma::mat prior_chol = arma::chol(prior_cov);
const arma::mat prior_rooti = arma::trans(arma::inv(arma::trimatu(prior_chol)));
const double prior_rootisum = arma::sum(log(prior_rooti.diag()));

//' @rdname mixtureexample_logprior
//' @title Evaluate mixture example prior density
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mixtureexample_logprior(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = -(static_cast<double>(dimension)/2.0) * log2pi + prior_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = prior_rooti * arma::trans(x.row(i) - prior_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname mixtureexample_gradlogprior
//' @title Evaluate gradient of mixture example prior density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat mixtureexample_gradlogprior(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (prior_mean - x.row(i)) * prior_precision;
  }
  
  return(output);
}

//' @rdname mixtureexample_partialderivative_logprior
//' @title Evaluate partial derivative of mixture example prior density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double mixtureexample_partialderivative_logprior(int dim, arma::rowvec x){
  arma::rowvec xc = prior_mean - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += prior_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname mixtureexample_sampleprior
//' @title Sample from mixture example prior distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat mixtureexample_sampleprior(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = prior_mean + Z * prior_chol;
  }
  
  return(output);
}

// Likelihood
arma::mat construct_like_obs(double xi){
  arma::mat output(ncomponents, dimension); 
  output(0, 0) = -xi;
  output(0, 1) = xi;
  
  output(1, 0) = xi;
  output(1, 1) = xi;
  
  output(2, 0) = -xi;
  output(2, 1) = -xi;
  
  output(3, 0) = xi;
  output(3, 1) = -xi;
  
  return(output);
  
}
const arma::mat like_obs = construct_like_obs(like_xi);

arma::mat construct_like_cov(double sign){
  arma::mat output = arma::zeros<arma::mat>(dimension, dimension);
  for (int i = 0; i < dimension; i++){
    for (int j = 0; j < dimension; j++){
      if (i == j) {
        output(i,j) = like_var;
      } else {
        output(i,j) = sign * like_rho;
      }
    }
  }
  return(output);
  
}
const arma::mat like_cov_minus = construct_like_cov(-1.0);
const arma::mat like_cov_plus = construct_like_cov(1.0);

const arma::mat like_precision_minus = arma::inv(like_cov_minus);
const arma::mat like_rooti_minus = arma::trans(arma::inv(arma::trimatu(arma::chol(like_cov_minus))));
const double like_rootisum_minus = arma::sum(log(like_rooti_minus.diag()));
const double constants_minus = -(static_cast<double>(dimension)/2.0) * log2pi + like_rootisum_minus;

const arma::mat like_precision_plus = arma::inv(like_cov_plus);
const arma::mat like_rooti_plus = arma::trans(arma::inv(arma::trimatu(arma::chol(like_cov_plus))));
const double like_rootisum_plus = arma::sum(log(like_rooti_plus.diag()));
const double constants_plus = -(static_cast<double>(dimension)/2.0) * log2pi + like_rootisum_plus;

//' @rdname mixtureexample_loglikelihood
//' @title Evaluate mixture example loglikelihood function
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mixtureexample_loglikelihood(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double maxlogdensities;
  arma::vec log_densities(ncomponents); 
  arma::vec densities(ncomponents);
  arma::vec z(dimension);
  
  for (int i = 0; i < n; i++){
    z = like_rooti_minus * arma::trans(x.row(i) - like_obs.row(0));
    log_densities(0) = std::log(like_weights(0)) + constants_minus - 0.5 * arma::sum(z % z);
    
    z = like_rooti_plus * arma::trans(x.row(i) - like_obs.row(1));
    log_densities(1) = std::log(like_weights(1)) + constants_plus - 0.5 * arma::sum(z % z);
    
    z = like_rooti_plus * arma::trans(x.row(i) - like_obs.row(2));
    log_densities(2) = std::log(like_weights(2)) + constants_plus - 0.5 * arma::sum(z % z);
    
    z = like_rooti_minus * arma::trans(x.row(i) - like_obs.row(3));
    log_densities(3) = std::log(like_weights(3)) + constants_minus - 0.5 * arma::sum(z % z);
    
    maxlogdensities = arma::max(log_densities);
    densities = arma::exp(log_densities - maxlogdensities);
    output(i) = std::log(sum(densities)) + maxlogdensities;
    
  }

  return(output);
}

//' @rdname mixtureexample_gradloglikelihood
//' @title Evaluate gradient of mixture example loglikelihood function
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat mixtureexample_gradloglikelihood(arma::mat x){
  int n = x.n_rows;
  double log_normaldensity, maxlogdensities, log_likelihood;
  arma::rowvec gradient(dimension);
  arma::mat output(n, dimension);
  arma::vec z(dimension);
  arma::vec log_densities(ncomponents); 
  arma::vec densities(ncomponents);

  for (int i = 0; i < n; i++){
    z = like_rooti_minus * arma::trans(x.row(i) - like_obs.row(0));
    log_normaldensity = constants_minus - 0.5 * arma::sum(z % z);
    log_densities(0) = std::log(like_weights(0)) + log_normaldensity;
    gradient = (like_weights(0) * std::exp(log_normaldensity)) * (like_obs.row(0) - x.row(i)) * like_precision_minus;
      
    z = like_rooti_plus * arma::trans(x.row(i) - like_obs.row(1));
    log_normaldensity = constants_plus - 0.5 * arma::sum(z % z);
    log_densities(1) = std::log(like_weights(1)) + log_normaldensity;
    gradient = gradient + (like_weights(1) * std::exp(log_normaldensity)) * (like_obs.row(1) - x.row(i)) * like_precision_plus;
    
    z = like_rooti_plus * arma::trans(x.row(i) - like_obs.row(2));
    log_normaldensity = constants_plus - 0.5 * arma::sum(z % z);
    log_densities(2) = std::log(like_weights(2)) + log_normaldensity;
    gradient = gradient + (like_weights(2) * std::exp(log_normaldensity)) * (like_obs.row(2) - x.row(i)) * like_precision_plus;
    
    z = like_rooti_minus * arma::trans(x.row(i) - like_obs.row(3));
    log_normaldensity = constants_minus - 0.5 * arma::sum(z % z);
    log_densities(3) = std::log(like_weights(3)) + log_normaldensity;
    gradient = gradient + (like_weights(3) * std::exp(log_normaldensity)) * (like_obs.row(3) - x.row(i)) * like_precision_minus;
    
    maxlogdensities = arma::max(log_densities);
    densities = arma::exp(log_densities - maxlogdensities);
    log_likelihood = std::log(sum(densities)) + maxlogdensities;
    gradient = gradient / std::exp(log_likelihood);
    output.row(i) = gradient;
    
  }

  return(output);
}

//' @rdname mixtureexample_partialderivative_loglikelihood
//' @title Evaluate partial derivative of mixture example likelihood function
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double mixtureexample_partialderivative_loglikelihood(int dim, arma::rowvec x){
  double output, log_normaldensity, partial_derivative, maxlogdensities, log_likelihood;
  arma::vec z(dimension);
  arma::vec log_densities(ncomponents); 
  arma::vec densities(ncomponents);
  arma::rowvec xc;
  
  z = like_rooti_minus * arma::trans(x - like_obs.row(0));
  log_normaldensity = constants_minus - 0.5 * arma::sum(z % z);
  log_densities(0) = std::log(like_weights(0)) + log_normaldensity;
  xc = like_obs.row(0) - x;
  partial_derivative = 0;
  for (int i = 0; i < dimension; i++){
    partial_derivative += like_precision_minus(dim, i) * xc(i);
  }
  output = like_weights(0) * std::exp(log_normaldensity) * partial_derivative;
  
  z = like_rooti_plus * arma::trans(x - like_obs.row(1));
  log_normaldensity = constants_plus - 0.5 * arma::sum(z % z);
  log_densities(1) = std::log(like_weights(1)) + log_normaldensity;
  xc = like_obs.row(1) - x;
  partial_derivative = 0;
  for (int i = 0; i < dimension; i++){
    partial_derivative += like_precision_plus(dim, i) * xc(i);
  }
  output += like_weights(1) * std::exp(log_normaldensity) * partial_derivative;
  
  z = like_rooti_plus * arma::trans(x - like_obs.row(2));
  log_normaldensity = constants_plus - 0.5 * arma::sum(z % z);
  log_densities(2) = std::log(like_weights(2)) + log_normaldensity;
  xc = like_obs.row(2) - x;
  partial_derivative = 0;
  for (int i = 0; i < dimension; i++){
    partial_derivative += like_precision_plus(dim, i) * xc(i);
  }
  output += like_weights(2) * std::exp(log_normaldensity) * partial_derivative;
  
  z = like_rooti_minus * arma::trans(x - like_obs.row(3));
  log_normaldensity = constants_minus - 0.5 * arma::sum(z % z);
  log_densities(3) = std::log(like_weights(3)) + log_normaldensity;
  xc = like_obs.row(3) - x;
  partial_derivative = 0;
  for (int i = 0; i < dimension; i++){
    partial_derivative += like_precision_minus(dim, i) * xc(i);
  }
  output += like_weights(3) * std::exp(log_normaldensity) * partial_derivative;
  
  maxlogdensities = arma::max(log_densities);
  densities = arma::exp(log_densities - maxlogdensities);
  log_likelihood = std::log(sum(densities)) + maxlogdensities;
  output /= std::exp(log_likelihood);

  return(output);
  
}

//' @rdname mixtureexample_posterior_cov
//' @title Compute posterior covariance of mixture example
//' @param component of mixture
//' @return posterior covariance
//' @export
// [[Rcpp::export]]
arma::mat mixtureexample_posterior_cov(int component){
  arma::mat output(dimension, dimension);
  
  switch (component) {
    case 0 : 
      output = arma::inv_sympd(prior_precision + like_precision_minus);
      break;
    case 1 : 
      output = arma::inv_sympd(prior_precision + like_precision_plus);
      break;
    case 2 : 
      output = arma::inv_sympd(prior_precision + like_precision_plus);
      break;
    case 3 : 
      output = arma::inv_sympd(prior_precision + like_precision_minus);
      break;
    
  }
  
  return(output);
}

//' @rdname mixtureexample_posterior_mean
//' @title Compute posterior mean of mixture example
//' @param component of mixture
//' @return posterior mean
//' @export
// [[Rcpp::export]]
arma::rowvec mixtureexample_posterior_mean(int component){
  arma::rowvec output(dimension);
  
  switch (component) {
    case 0 : 
      output = (prior_mean * prior_precision + like_obs.row(0) * like_precision_minus) * mixtureexample_posterior_cov(0);
      break;
    case 1 : 
      output = (prior_mean * prior_precision + like_obs.row(1) * like_precision_plus) * mixtureexample_posterior_cov(1);
      break;  
    case 2 : 
      output = (prior_mean * prior_precision + like_obs.row(2) * like_precision_plus) * mixtureexample_posterior_cov(2);
      break;
    case 3 : 
      output = (prior_mean * prior_precision + like_obs.row(3) * like_precision_minus) * mixtureexample_posterior_cov(3);
      break;
    
  }
  
  return(output);
}


//' @rdname mixtureexample_log_normconst
//' @title Compute log normalizing constant of mixture example
//' @return log normalizing constant
//' @export
// [[Rcpp::export]]
double mixtureexample_log_normconst(){
  arma::rowvec posterior_mean = mixtureexample_posterior_mean(0);
  arma::mat posterior_cov = mixtureexample_posterior_cov(0);
  arma::mat posterior_precision = prior_precision + like_precision_minus;
  double logdet_priorcov = log(arma::det(prior_cov));
  double logdet_likecov = log(arma::det(like_cov_minus));
  double logdet_posteriorcov = log(arma::det(posterior_cov));
  double output = - (static_cast<double>(dimension)/2.0) * log2pi - 0.5 * logdet_likecov -
    0.5 * logdet_priorcov + 0.5 * logdet_posteriorcov -
    0.5 * arma::sum( prior_mean % (prior_mean * prior_precision) ) -
    0.5 * arma::sum( like_obs.row(0) %  (like_obs.row(0) * like_precision_minus) ) +
    0.5 * sum( posterior_mean % (posterior_mean * posterior_precision) );

  return(output);
}

//' @rdname mixtureexample_logposterior
//' @title Evaluate mixture example posterior density
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mixtureexample_logposterior(arma::mat x){
  
  arma::mat posterior_means(ncomponents, dimension);
  for (int i = 0; i < ncomponents; i++){
    posterior_means.row(i) = mixtureexample_posterior_mean(i); 
  }
  
  arma::mat posterior_precision_minus = prior_precision + like_precision_minus;
  arma::mat posterior_cov_minus = arma::inv_sympd(posterior_precision_minus);
  arma::mat posterior_chol_minus = arma::chol(posterior_cov_minus);
  arma::mat posterior_rooti_minus = arma::trans(arma::inv(arma::trimatu(posterior_chol_minus)));
  double posterior_rootisum_minus = arma::sum(log(posterior_rooti_minus.diag()));
  double posterior_constants_minus = -(static_cast<double>(dimension)/2.0) * log2pi + posterior_rootisum_minus;
  
  arma::mat posterior_precision_plus = prior_precision + like_precision_plus;
  arma::mat posterior_cov_plus = arma::inv_sympd(posterior_precision_plus);
  arma::mat posterior_chol_plus = arma::chol(posterior_cov_plus);
  arma::mat posterior_rooti_plus = arma::trans(arma::inv(arma::trimatu(posterior_chol_plus)));
  double posterior_rootisum_plus = arma::sum(log(posterior_rooti_plus.diag()));
  double posterior_constants_plus = -(static_cast<double>(dimension)/2.0) * log2pi + posterior_rootisum_plus;
                                               
  int n = x.n_rows;
  arma::vec output(n);
  double maxlogdensities;
  arma::vec log_densities(ncomponents); 
  arma::vec densities(ncomponents);
  arma::vec z(dimension);
  
  for (int i = 0; i < n; i++){
    
    z = posterior_rooti_minus * arma::trans(x.row(i) - posterior_means.row(0));
    log_densities(0) = std::log(like_weights(0)) + posterior_constants_minus - 0.5 * arma::sum(z % z);
    
    z = posterior_rooti_plus * arma::trans(x.row(i) - posterior_means.row(1));
    log_densities(1) = std::log(like_weights(1)) + posterior_constants_plus - 0.5 * arma::sum(z % z);
    
    z = posterior_rooti_plus * arma::trans(x.row(i) - posterior_means.row(2));
    log_densities(2) = std::log(like_weights(2)) + posterior_constants_plus - 0.5 * arma::sum(z % z);
    
    z = posterior_rooti_minus * arma::trans(x.row(i) - posterior_means.row(3));
    log_densities(3) = std::log(like_weights(3)) + posterior_constants_minus - 0.5 * arma::sum(z % z);
    
    maxlogdensities = arma::max(log_densities);
    densities = arma::exp(log_densities - maxlogdensities);
    output(i) = std::log(sum(densities)) + maxlogdensities;
    
  }

  return(output);
}

// Gibbs flow computation
const double lowerbound = 10.0;
const double upperbound = -10.0;
const int ngridpoints = 200;
const arma::vec grid = arma::linspace(lowerbound, upperbound, ngridpoints);
const double mesh = grid(1) - grid(0);

//' @rdname mixtureexample_gibbsflow
//' @title Compute Gibbs flow for mixture example
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
Rcpp::List mixtureexample_gibbsflow(double stepsize, double lambda, double derivative_lambda,
                              arma::mat xparticles, arma::vec logdensity){
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
      loglikelihood_x = as_scalar(mixtureexample_loglikelihood(xpoint));
      target_x = std::exp( as_scalar(mixtureexample_logprior(xpoint)) + lambda * loglikelihood_x );
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
        loglikelihood_gridpoint = as_scalar(mixtureexample_loglikelihood(gridpoint));
        target_gridpoint = std::exp( as_scalar(mixtureexample_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
        loglikelihood_gridpoint = as_scalar(mixtureexample_loglikelihood(gridpoint));
        target_gridpoint = std::exp( as_scalar(mixtureexample_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
      logtarget_partial_derivative = mixtureexample_partialderivative_logprior(i, xpoint) + lambda * mixtureexample_partialderivative_loglikelihood(i, xpoint);
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

//' @rdname mixtureexample_gibbsvelocity
//' @title Compute Gibbs flow for mixture example
//' @param time time
//' @param xpoint particle position
//' @param exponent exponent of tempering schedule
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::rowvec mixtureexample_gibbsvelocity(double time, arma::rowvec xpoint, double exponent){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * pow(time, exponent - 1.0);
  
  // Output 
  arma::rowvec gibbs_velocity(dimension);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double loglikelihood_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  arma::rowvec gridpoint;
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
    
  loglikelihood_x = as_scalar(mixtureexample_loglikelihood(xpoint));
  target_x = std::exp( as_scalar(mixtureexample_logprior(xpoint)) + lambda * loglikelihood_x );
  loglikelihood_target_x = loglikelihood_x * target_x;
  for (int i = 0; i < dimension; i++){
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
      loglikelihood_gridpoint = as_scalar(mixtureexample_loglikelihood(gridpoint));
      target_gridpoint = std::exp( as_scalar(mixtureexample_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
      loglikelihood_gridpoint = as_scalar(mixtureexample_loglikelihood(gridpoint));
      target_gridpoint = std::exp( as_scalar(mixtureexample_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
    gibbs_velocity(i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
      
  }
    
  return(gibbs_velocity);
    
}





