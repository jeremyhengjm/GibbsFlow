#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Problem settings
const int ngrid = 20;
const int dimension = ngrid * ngrid;
double double_dimension = static_cast<double>(dimension);
const double parameter_sigmasq = 1.91;
const double parameter_mu = log(126.0) - 0.5 * parameter_sigmasq;
const double parameter_inverse_beta = 33;
const double parameter_area = 1.0 / double_dimension;

// Prior 
const double log2pi = std::log(2.0 * arma::datum::pi);
const arma::rowvec prior_mean = parameter_mu * arma::ones<arma::rowvec>(dimension);
arma::mat construct_prior_cov(int dimension){
  arma::mat prior_cov = arma::zeros<arma::mat>(dimension, dimension);
  double double_ngrid = static_cast<double>(ngrid);
  double double_m, double_n;
  double index_m1, index_m2; 
  double index_n1, index_n2;
  double distance;
  for (int m = 0; m < dimension; m++){
    for (int n = 0; n < dimension; n++){
      double_m = static_cast<double>(m);
      double_n = static_cast<double>(n);
      
      index_m1 = std::floor(double_m / double_ngrid) + 1;
      index_n1 = std::floor(double_n / double_ngrid) + 1;
      
      index_m2 = std::fmod(double_m, double_ngrid) + 1;
      index_n2 = std::fmod(double_n, double_ngrid) + 1;  
      
      distance = sqrt( std::pow(index_m1 - index_n1, 2.0) + std::pow(index_m2 - index_n2, 2.0) );
      prior_cov(m, n) = parameter_sigmasq * exp(- parameter_inverse_beta * distance / double_ngrid);
      
    }
  }
  return(prior_cov);
}
const arma::mat prior_cov = construct_prior_cov(dimension);
const arma::mat prior_precision = arma::inv_sympd(prior_cov);
const arma::mat prior_chol = arma::chol(prior_cov);
const arma::mat prior_rooti = arma::trans(arma::inv(arma::trimatu(prior_chol)));
const double prior_rootisum = arma::sum(log(prior_rooti.diag()));

//' @rdname coxprocess_logprior
//' @title Evaluate multivariate Gaussian prior density
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_logprior(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = - 0.5 * double_dimension * log2pi + prior_rootisum;

  for (int i = 0; i < n; i++){
    arma::vec z = prior_rooti * arma::trans(x.row(i) - prior_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }

  return(output);
}

//' @rdname coxprocess_gradlogprior
//' @title Evaluate gradient of multivariate Gaussian prior density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gradlogprior(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (prior_mean - x.row(i)) * prior_precision;
  }
  
  return(output);
}

//' @rdname coxprocess_partialderivative_logprior
//' @title Evaluate partial derivative of multivariate Gaussian prior density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point  
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double coxprocess_partialderivative_logprior(int dim, arma::rowvec x){
  arma::rowvec xc = prior_mean - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += prior_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname coxprocess_sampleprior
//' @title Sample from multivariate Gaussian prior distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_sampleprior(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = prior_mean + Z * prior_chol;
  }
  
  return(output);
}

arma::mat precompute_conditional_term(int dimension){
  arma::mat output = arma::zeros<arma::mat>(dimension, dimension-1);
  arma::mat inv_submat(dimension-1, dimension-1);
  arma::vec index = arma::regspace(0, dimension-1);
  for (int i = 0; i < dimension; i++){
    arma::uvec cond_index = arma::find(index == i);
    arma::uvec other_index = arma::find(index != i);
    inv_submat = arma::inv_sympd(prior_cov.submat(other_index, other_index));
    output.row(i) = prior_cov.submat(cond_index, other_index) * inv_submat;
  }
  return(output);
}
const arma::mat conditional_term = precompute_conditional_term(dimension);

arma::vec precompute_conditional_var(int dimension){
  arma::vec output(dimension);
  arma::vec index = arma::regspace(0, dimension-1);
  for (int i = 0; i < dimension; i++){
    arma::uvec cond_index = arma::find(index == i);
    arma::uvec other_index = arma::find(index != i);
    output(i) = prior_cov(i, i) - as_scalar(conditional_term.row(i) * prior_cov.submat(other_index, cond_index));
  }
  return(output);
}
const arma::vec conditional_var = precompute_conditional_var(dimension);

//' @rdname coxprocess_logprior_conditional
//' @title Evaluate conditional density of multivariate Gaussian prior 
//' @param dim dimension of interest (dim = 0 is the first coordinate)
//' @param x evaluation points along dim
//' @param other_x conditional values
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_logprior_conditional(int dim, arma::vec x, arma::vec other_x){
  int n = x.n_elem;
  arma::vec output(n);
  arma::vec conditional_centered_mean = other_x - parameter_mu;
  double conditional_mean = parameter_mu + as_scalar(conditional_term.row(dim) * conditional_centered_mean);
  double constants = - 0.5 * ( log2pi + log(conditional_var(dim)) );
  double conditional_precision = 1.0 / conditional_var(dim);
  double centered_mean;
  
  for (int i = 0; i < n; i++){
    centered_mean = x(i) - conditional_mean;
    output(i) = constants - 0.5 * conditional_precision * centered_mean * centered_mean;
  }
  
  return(output);
}

// Likelihood

//' @rdname coxprocess_loglikelihood
//' @title Evaluate Cox process loglikelihood function
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_loglikelihood(arma::mat x, arma::vec counts){
  int N = x.n_rows;
  arma::vec output(N);
  double cumsum;
  
  for (int n = 0; n < N; n++){
    cumsum = 0;
    for (int i = 0; i < dimension; i++){
      cumsum += x(n, i) * counts(i) - parameter_area * exp(x(n, i));
    }
    output(n) = cumsum;
  }
  
  return(output);
}

//' @rdname coxprocess_loglikelihood_term
//' @title Evaluate a particular term of Cox process loglikelihood function
//' @param dim loglikelihood term depends only on this dimension (dim = 0 is the first coordinate)
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_loglikelihood_term(int dim, arma::vec x, arma::vec counts){
  int N = x.n_rows;
  arma::vec output(N);
  
  for (int n = 0; n < N; n++){
    output(n) = x(n) * counts(dim) - parameter_area * exp(x(n));
  }
  
  return(output);
}

//' @rdname coxprocess_gradloglikelihood
//' @title Evaluate gradient of Cox process loglikelihood function
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gradloglikelihood(arma::mat x, arma::vec counts){
  int N = x.n_rows;
  arma::mat output(N, dimension);
  
  for (int n = 0; n < N; n++){
    for (int i = 0; i < dimension; i++){
      output(n, i) = counts(i) - parameter_area * exp(x(n, i));
    }
  }
  
  return(output);
}

//' @rdname coxprocess_partialderivative_loglikelihood
//' @title Evaluate partial derivative of Cox process likelihood function
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point  
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double coxprocess_partialderivative_loglikelihood(int dim, arma::rowvec x, arma::vec counts){
  double output = counts(dim) - parameter_area * exp(x(dim));

  return(output);
}

// Gibbs flow computation
// const double lowerbound = -4.0;
// const double upperbound = 12.0;
const double lowerbound = parameter_mu - 6.0 * std::sqrt(parameter_sigmasq);
const double upperbound = parameter_mu + 6.0 * std::sqrt(parameter_sigmasq);
const int ngridpoints = 40;
const arma::vec grid = arma::linspace(lowerbound, upperbound, ngridpoints);
const double mesh = grid(1) - grid(0);

//' @rdname coxprocess_gibbsflow
//' @title Compute Gibbs flow for Cox process model
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @param counts dataset 
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List coxprocess_gibbsflow(double stepsize, double lambda, double derivative_lambda,
                                    arma::mat xparticles, arma::vec logdensity, arma::vec counts){
  int N = xparticles.n_rows; // no. of particles
  
  // Output
  arma::rowvec log_jacobian_dets(N);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double log_jacobian_det;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double target_x, loglikelihood_target_x, loglikelihood_gridpoint_i, loglikelihood_gridpoint;  
  double target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double logprior_conditional_x, logprior_conditional_gridpoint;
  
  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean;
  
  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    log_jacobian_det = 0;
    
    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }
    
    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;
      
      // Evaluate conditional density 
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      target_x = std::exp(logprior_conditional_x + lambda * loglikelihood_x_i); 
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
      
      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional density 
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        target_gridpoint = std::exp(logprior_conditional_gridpoint + lambda * loglikelihood_gridpoint_i);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = loglikelihood_gridpoint * target_gridpoint;
      }
      
      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional density 
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        target_gridpoint = std::exp(logprior_conditional_gridpoint + lambda * loglikelihood_gridpoint_i);
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
      logtarget_partial_derivative = coxprocess_partialderivative_logprior(i, xpoint) + lambda * coxprocess_partialderivative_loglikelihood(i, xpoint, counts);
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - loglikelihood_x) - 
        gibbs_velocity * logtarget_partial_derivative;
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
      xpoint(i) = xpoint(i) + stepsize * gibbs_velocity;
      
      // Update loglikelihood
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x = loglikelihood_x_i + loglikelihood_x_lessi;
      
    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;
    
  }
  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles, 
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);
  
}

//' @rdname coxprocess_gibbsvelocity
//' @title Compute Gibbs velocity for Cox process model
//' @param time time
//' @param xparticles particle positions
//' @param exponent exponent of tempering schedule
//' @param counts dataset 
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gibbsvelocity(double time, arma::mat xparticles, double exponent, arma::vec counts){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * std::pow(time, exponent - 1.0);
  int N = xparticles.n_rows; // no. of particles
  
  // Output 
  arma::mat gibbs_velocity(N, dimension);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double target_x, loglikelihood_target_x, loglikelihood_gridpoint_i, loglikelihood_gridpoint;  
  double target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double logprior_conditional_x, logprior_conditional_gridpoint;
  
  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean;
  
  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    
    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }
    
    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;
      
      // Evaluate conditional density 
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      target_x = std::exp(logprior_conditional_x + lambda * loglikelihood_x_i); 
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
      
      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional density 
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        target_gridpoint = std::exp(logprior_conditional_gridpoint + lambda * loglikelihood_gridpoint_i);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = loglikelihood_gridpoint * target_gridpoint;
      }
      
      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional density 
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        target_gridpoint = std::exp(logprior_conditional_gridpoint + lambda * loglikelihood_gridpoint_i);
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
      gibbs_velocity(n, i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
      
    }
  }
  
  return(gibbs_velocity);  
  
}

//' @rdname coxprocess_stein_variational_importance_sampling
//' @title Perform Stein variational importance sampling
//' @param nfparticles number of follower particles
//' @param nlparticles number of leader particles
//' @param niterations number of time iterations 
//' @param stepsize stepsize of each iteration
//' @param counts dataset
//' @return list with keys:
//' \code{xparticles} new particle locations
//' @export
// [[Rcpp::export]]
Rcpp::List coxprocess_stein_variational_importance_sampling(int nfparticles, int nlparticles, int niterations, double stepsize, arma::vec counts){
  int counter;
  double double_nlparticles = static_cast<double>(nlparticles);
  double difference, distance, bandwidth;
  double drift, partial_logtarget, kernel_value, grad_kernel, log_jacobian_det, log_jacobian_det_temp;
  double maxlogweights, ess, log_normconst;
  arma::mat dist_xl_xl(nlparticles, nlparticles);
  arma::vec dist_xl_xl_vec(nlparticles * nlparticles);
  arma::mat dist_xl_xf(nlparticles, nfparticles);
  arma::mat xlparticles = coxprocess_sampleprior(nlparticles);
  arma::mat xfparticles = coxprocess_sampleprior(nfparticles);
  arma::vec logproposal = coxprocess_logprior(xfparticles);
  arma::vec logweights(nfparticles);
  arma::vec weights(nfparticles);
  arma::vec normweights(nfparticles);
  arma::vec ess_percentage(niterations);

  for (int l = 0; l < niterations; l++){
    Rcpp::Rcout << "Iterations:" << std::endl << l << std::endl;
    // Compute pairwise distance
    counter = 0;
    for (int j = 0; j < nlparticles; j++){
      // Between leader particles
      for (int p = 0; p < nlparticles; p++){
        distance = 0;
        for (int i = 0; i < dimension; i++){
          difference = xlparticles(j, i) - xlparticles(p, i);
          distance += difference * difference;
        }
        dist_xl_xl(j, p) = distance;
        dist_xl_xl_vec(counter) = distance;
        counter += 1;
      }
      
      // // Between leader and follower particles
      for (int n = 0; n < nfparticles; n++){
        distance = 0;
        for (int i = 0; i < dimension; i++){
          difference = xlparticles(j, i) - xfparticles(n, i);
          distance += difference * difference;
        }
        dist_xl_xf(j, n) = distance;
      }
      
    }
    
    // Set bandwidth
    bandwidth = 0.5 * std::pow(arma::median(dist_xl_xl_vec), 2) / std::log(double_nlparticles + 1.0);
    
    // Move follower particles
    for (int n = 0; n < nfparticles; n++){
      log_jacobian_det = 0;

      for (int i = 0; i < dimension; i++){
        drift = 0;
        log_jacobian_det_temp = 0;
        for (int j = 0; j < nlparticles; j++){
          partial_logtarget = coxprocess_partialderivative_logprior(i, xlparticles.row(j)) +
                              coxprocess_partialderivative_loglikelihood(i, xlparticles.row(j), counts);
          kernel_value = std::exp(- dist_xl_xf(j, n) / bandwidth);
          grad_kernel = 2.0 * (xlparticles(j, i) - xfparticles(n, i)) * kernel_value / bandwidth;
          drift += partial_logtarget * kernel_value + grad_kernel;
          log_jacobian_det_temp += partial_logtarget * grad_kernel + grad_kernel * grad_kernel / kernel_value - 2.0 * kernel_value / bandwidth;
        }
        // Rcpp::Rcout << "Drift:" << std::endl << drift << std::endl;
        xfparticles(n, i) = xfparticles(n, i) + stepsize * drift / double_nlparticles;
        log_jacobian_det -= std::log(std::abs(1.0 + stepsize * log_jacobian_det_temp / double_nlparticles));

      }
      logproposal(n) += log_jacobian_det;
      logweights(n) = as_scalar(coxprocess_logprior(xfparticles.row(n))) + as_scalar(coxprocess_loglikelihood(xfparticles.row(n), counts)) - logproposal(n);
    }
    
    // // Weight
    maxlogweights = arma::max(logweights);
    weights = arma::exp(logweights - maxlogweights);
    normweights = weights / arma::sum(weights);
    ess = 1.0 / arma::sum(arma::square(normweights));
    ess_percentage(l) = ess / nfparticles * 100;
    Rcpp::Rcout << "ESS%:" << std::endl << ess_percentage(l) << std::endl;
    log_normconst = std::log(arma::mean(weights)) + maxlogweights;
    
    // Move leader particles 
    for (int n = 0; n < nlparticles; n++){
      for (int i = 0; i < dimension; i++){
        drift = 0;
        for (int j = 0; j < nlparticles; j++){
          partial_logtarget = coxprocess_partialderivative_logprior(i, xlparticles.row(j)) +
                              coxprocess_partialderivative_loglikelihood(i, xlparticles.row(j), counts);
          kernel_value = std::exp(- dist_xl_xl(j, n) / bandwidth);
          grad_kernel = 2.0 * (xlparticles(j, i) - xlparticles(n, i)) * kernel_value / bandwidth;
          drift += partial_logtarget * kernel_value + grad_kernel;
        }
        xlparticles(n, i) = xlparticles(n, i) + stepsize * drift / double_nlparticles;
      }
    }
    

  }

  return Rcpp::List::create(Rcpp::Named("xparticles") = xfparticles,
                            Rcpp::Named("ess") = ess_percentage,
                            Rcpp::Named("log_normconst") = log_normconst);
  
}

// Load EP Gaussian approximation
arma::rowvec load_ep_mean(){
  arma::rowvec output;
  output.load("inst/coxprocess/ep_mean.bin");
  return(output);
}
const arma::rowvec ep_mean = load_ep_mean();

arma::mat load_ep_cov(){
  arma::mat output;
  output.load("inst/coxprocess/ep_cov.bin");
  return(output);
}
const arma::mat ep_cov = load_ep_cov();
// const arma::mat proposal_cov = load_ep_cov();
// const arma::mat ep_cov = arma::diagmat(proposal_cov);
const arma::mat ep_precision  = arma::inv_sympd(ep_cov);
const arma::mat ep_chol = arma::chol(ep_cov);
const arma::mat ep_rooti = arma::trans(arma::inv(arma::trimatu(ep_chol)));
const double ep_rootisum = arma::sum(log(ep_rooti.diag()));

//' @rdname coxprocess_log_ep_proposal
//' @title Evaluate EP Gaussian proposal density
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_log_ep_proposal(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = - 0.5 * double_dimension * log2pi + ep_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = ep_rooti * arma::trans(x.row(i) - ep_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname coxprocess_gradlog_ep_proposal
//' @title Evaluate gradient of EP Gaussian proposal density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gradlog_ep_proposal(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (ep_mean - x.row(i)) * ep_precision;
  }
  
  return(output);
}

//' @rdname coxprocess_sample_ep_proposal
//' @title Sample from EP Gaussian proposal distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_sample_ep_proposal(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = ep_mean + Z * ep_chol;
  }
  
  return(output);
}

arma::mat precompute_ep_conditional_term(int dimension){
  arma::mat output = arma::zeros<arma::mat>(dimension, dimension-1);
  arma::mat inv_submat(dimension-1, dimension-1);
  arma::vec index = arma::regspace(0, dimension-1);
  for (int i = 0; i < dimension; i++){
    arma::uvec cond_index = arma::find(index == i);
    arma::uvec other_index = arma::find(index != i);
    inv_submat = arma::inv_sympd(ep_cov.submat(other_index, other_index));
    output.row(i) = ep_cov.submat(cond_index, other_index) * inv_submat;
  }
  return(output);
}

const arma::mat ep_conditional_term = precompute_ep_conditional_term(dimension);

arma::mat precompute_ep_conditional_var(int dimension){
  arma::vec output(dimension);
  arma::vec index = arma::regspace(0, dimension-1);
  for (int i = 0; i < dimension; i++){
    arma::uvec cond_index = arma::find(index == i);
    arma::uvec other_index = arma::find(index != i);
    output(i) = ep_cov(i, i) - as_scalar(ep_conditional_term.row(i) * ep_cov.submat(other_index, cond_index));
  }
  return(output);
}

const arma::vec ep_conditional_var = precompute_ep_conditional_var(dimension);

//' @rdname coxprocess_log_ep_proposal_conditional
//' @title Evaluate conditional density of EP Gaussian proposal 
//' @param dim dimension of interest (dim = 0 is the first coordinate)
//' @param x evaluation points along dim
//' @param other_x conditional values
//' @param other_index logical vector indexing conditioned coordinates
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_log_ep_proposal_conditional(int dim, arma::vec x, arma::vec other_x, arma::uvec other_index){
  int n = x.n_elem;
  arma::vec output(n);
  
  arma::vec conditional_centered_mean = other_x - ep_mean.elem(other_index);
  
  double conditional_mean = ep_mean(dim) + as_scalar(ep_conditional_term.row(dim) * conditional_centered_mean);
  double constants = - 0.5 * ( log2pi + log(ep_conditional_var(dim)) );
  double conditional_precision = 1.0 / ep_conditional_var(dim);
  double centered_mean;

  for (int i = 0; i < n; i++){
    centered_mean = x(i) - conditional_mean;
    output(i) = constants - 0.5 * conditional_precision * centered_mean * centered_mean;
  }

  return(output);
}

//' @rdname coxprocess_partialderivative_log_ep_proposal
//' @title Evaluate partial derivative of EP Gaussian proposal density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point  
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double coxprocess_partialderivative_log_ep_proposal(int dim, arma::rowvec x){ 
  arma::rowvec xc = ep_mean - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += ep_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname coxprocess_gibbsflow_ep_proposal
//' @title Compute Gibbs flow for Cox process model with EP Gaussian proposal initialization
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @param counts dataset
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List coxprocess_gibbsflow_ep_proposal(double stepsize, double lambda, double derivative_lambda,
                                            arma::mat xparticles, arma::vec logdensity, arma::vec counts){

  int N = xparticles.n_rows; // no. of particles

  // Output
  arma::rowvec log_jacobian_dets(N);

  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double log_jacobian_det;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double logprior_conditional_x, logproposal_conditional_x, log_artificial_likelihood_x;
  double logtarget_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, loglikelihood_gridpoint_i;
  double logprior_conditional_gridpoint, logproposal_conditional_gridpoint, log_artificial_likelihood_gridpoint;
  double logtarget_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double ep_conditional_mean, ep_conditional_constants, ep_conditional_precision, ep_centered_mean;

  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean, ep_conditional_centered_mean;

  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    log_jacobian_det = 0;

    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }

    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Precompute EP conditional distribution 
      ep_conditional_centered_mean = xpoint_lessi - ep_mean.elem(other_index);
      ep_conditional_mean = ep_mean(i) + as_scalar(ep_conditional_term.row(i) * ep_conditional_centered_mean);
      ep_conditional_constants = - 0.5 * ( log2pi + log(ep_conditional_var(i)) );
      ep_conditional_precision = 1.0 / ep_conditional_var(i);
      
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;

      // Evaluate conditional densities
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      ep_centered_mean = xpoint(i) - ep_conditional_mean;
      logproposal_conditional_x = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
      log_artificial_likelihood_x = logprior_conditional_x + loglikelihood_x - logproposal_conditional_x;
      logtarget_x = (1.0 - lambda) * logproposal_conditional_x + lambda * (logprior_conditional_x + loglikelihood_x_i);
      target_x = std::exp(logtarget_x);
      loglikelihood_target_x = log_artificial_likelihood_x * target_x;

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

      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);

      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;

      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);

        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;

        // Evaluate conditional densities
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        ep_centered_mean = lowergrid(p) - ep_conditional_mean;
        logproposal_conditional_gridpoint = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
      }

      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);

      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);

        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;

        // Evaluate conditional densities
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        ep_centered_mean = uppergrid(p) - ep_conditional_mean;
        logproposal_conditional_gridpoint = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_uppergrid(p) = target_gridpoint;
        loglikelihood_target_uppergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
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
      logtarget_partial_derivative = (1.0 - lambda) * coxprocess_partialderivative_log_ep_proposal(i, xpoint) + 
        lambda * ( coxprocess_partialderivative_logprior(i, xpoint) + coxprocess_partialderivative_loglikelihood(i, xpoint, counts) );
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - log_artificial_likelihood_x) -
        gibbs_velocity * logtarget_partial_derivative;
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
      xpoint(i) = xpoint(i) + stepsize * gibbs_velocity;

      // Update loglikelihood
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x = loglikelihood_x_i + loglikelihood_x_lessi;

    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;

  }
  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles,
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);

}

//' @rdname coxprocess_gibbsvelocity_ep_proposal
//' @title Compute Gibbs velocity for Cox process model with EP Gaussian proposal initialization
//' @param time time
//' @param xparticles particle positions
//' @param exponent exponent of tempering schedule
//' @param counts dataset
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gibbsvelocity_ep_proposal(double time, arma::mat xparticles, double exponent, arma::vec counts){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * std::pow(time, exponent - 1.0);
  int N = xparticles.n_rows; // no. of particles
  
  // Output
  arma::mat gibbs_velocity(N, dimension);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double logprior_conditional_x, logproposal_conditional_x, log_artificial_likelihood_x;
  double logtarget_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, loglikelihood_gridpoint_i;
  double logprior_conditional_gridpoint, logproposal_conditional_gridpoint, log_artificial_likelihood_gridpoint;
  double logtarget_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double ep_conditional_mean, ep_conditional_constants, ep_conditional_precision, ep_centered_mean;
  
  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean, ep_conditional_centered_mean;
  
  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    
    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }
    
    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Precompute EP conditional distribution 
      ep_conditional_centered_mean = xpoint_lessi - ep_mean.elem(other_index);
      ep_conditional_mean = ep_mean(i) + as_scalar(ep_conditional_term.row(i) * ep_conditional_centered_mean);
      ep_conditional_constants = - 0.5 * ( log2pi + log(ep_conditional_var(i)) );
      ep_conditional_precision = 1.0 / ep_conditional_var(i);
      
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;
      
      // Evaluate conditional densities
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      ep_centered_mean = xpoint(i) - ep_conditional_mean;
      logproposal_conditional_x = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
      log_artificial_likelihood_x = logprior_conditional_x + loglikelihood_x - logproposal_conditional_x;
      logtarget_x = (1.0 - lambda) * logproposal_conditional_x + lambda * (logprior_conditional_x + loglikelihood_x_i);
      target_x = std::exp(logtarget_x);
      loglikelihood_target_x = log_artificial_likelihood_x * target_x;
      
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
      
      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        ep_centered_mean = lowergrid(p) - ep_conditional_mean;
        logproposal_conditional_gridpoint = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
      }
      
      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        ep_centered_mean = uppergrid(p) - ep_conditional_mean;
        logproposal_conditional_gridpoint = ep_conditional_constants - 0.5 * ep_conditional_precision * ep_centered_mean * ep_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_uppergrid(p) = target_gridpoint;
        loglikelihood_target_uppergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
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
      gibbs_velocity(n, i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
      
    }
    
  }
  
  return(gibbs_velocity);
  
}

// Load VI Gaussian approximation
arma::rowvec load_vi_mean(){
  arma::rowvec output;
  output.load("inst/coxprocess/vi_mean.bin");
  return(output);
}
const arma::rowvec vi_mean = load_vi_mean();

arma::rowvec load_vi_var(){
  arma::rowvec output;
  output.load("inst/coxprocess/vi_var.bin");
  return(output);
}
const arma::rowvec vi_var = load_vi_var();
const arma::rowvec vi_std = arma::sqrt(vi_var);

//' @rdname coxprocess_log_vi_proposal
//' @title Evaluate VI Gaussian proposal density
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec coxprocess_log_vi_proposal(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = - 0.5 * ( double_dimension * log2pi + arma::sum(arma::log(vi_var)) );
  
  for (int i = 0; i < n; i++){
    output(i) = constants - 0.5 * arma::sum(arma::square(x.row(i) - vi_mean) / vi_var);
  }
  
  return(output);
}

//' @rdname coxprocess_gradlog_vi_proposal
//' @title Evaluate gradient of VI Gaussian proposal density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gradlog_vi_proposal(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (vi_mean - x.row(i)) / vi_var;
  }
  
  return(output);
}

//' @rdname coxprocess_sample_vi_proposal
//' @title Sample from VI Gaussian proposal distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_sample_vi_proposal(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = vi_mean + Z % vi_std;
  }
  
  return(output);
}

//' @rdname coxprocess_partialderivative_log_vi_proposal
//' @title Evaluate partial derivative of VI Gaussian proposal density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point  
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double coxprocess_partialderivative_log_vi_proposal(int dim, arma::rowvec x){ 
  double output = (vi_mean(dim) - x(dim)) / vi_var(dim);
  return(output);
  
}

//' @rdname coxprocess_gibbsflow_vi_proposal
//' @title Compute Gibbs flow for Cox process model with VI Gaussian proposal initialization
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @param counts dataset
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List coxprocess_gibbsflow_vi_proposal(double stepsize, double lambda, double derivative_lambda,
                                            arma::mat xparticles, arma::vec logdensity, arma::vec counts){
  
  int N = xparticles.n_rows; // no. of particles
  
  // Output
  arma::rowvec log_jacobian_dets(N);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double log_jacobian_det;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double logprior_conditional_x, logproposal_conditional_x, log_artificial_likelihood_x;
  double logtarget_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, loglikelihood_gridpoint_i;
  double logprior_conditional_gridpoint, logproposal_conditional_gridpoint, log_artificial_likelihood_gridpoint;
  double logtarget_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double vi_conditional_constants, vi_conditional_precision, vi_centered_mean;
  
  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean, vi_conditional_centered_mean;
  
  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    log_jacobian_det = 0;
    
    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }
    
    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Precompute VI conditional distribution 
      vi_conditional_constants = - 0.5 * ( log2pi + log(vi_var(i)) );
      vi_conditional_precision = 1.0 / vi_var(i);
        
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;
      
      // Evaluate conditional densities
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      vi_centered_mean = xpoint(i) - vi_mean(i);
      logproposal_conditional_x = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
      log_artificial_likelihood_x = logprior_conditional_x + loglikelihood_x - logproposal_conditional_x;
      logtarget_x = (1.0 - lambda) * logproposal_conditional_x + lambda * (logprior_conditional_x + loglikelihood_x_i);
      target_x = std::exp(logtarget_x);
      loglikelihood_target_x = log_artificial_likelihood_x * target_x;
      
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
      
      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        vi_centered_mean = lowergrid(p) - vi_mean(i);
        logproposal_conditional_gridpoint = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
      }
      
      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        vi_centered_mean = uppergrid(p) - vi_mean(i);
        logproposal_conditional_gridpoint = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_uppergrid(p) = target_gridpoint;
        loglikelihood_target_uppergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
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
      logtarget_partial_derivative = (1.0 - lambda) * coxprocess_partialderivative_log_vi_proposal(i, xpoint) + 
        lambda * ( coxprocess_partialderivative_logprior(i, xpoint) + coxprocess_partialderivative_loglikelihood(i, xpoint, counts) );
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - log_artificial_likelihood_x) -
        gibbs_velocity * logtarget_partial_derivative;
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
      xpoint(i) = xpoint(i) + stepsize * gibbs_velocity;
      
      // Update loglikelihood
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x = loglikelihood_x_i + loglikelihood_x_lessi;
      
    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;
    
  }
  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles,
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);
  
}

//' @rdname coxprocess_gibbsvelocity_vi_proposal
//' @title Compute Gibbs velocity for Cox process model with VI Gaussian proposal initialization
//' @param time time
//' @param xparticles particle positions
//' @param exponent exponent of tempering schedule
//' @param counts dataset
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::mat coxprocess_gibbsvelocity_vi_proposal(double time, arma::mat xparticles, double exponent, arma::vec counts){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * std::pow(time, exponent - 1.0);
  int N = xparticles.n_rows; // no. of particles
  
  // Output
  arma::mat gibbs_velocity(N, dimension);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double loglikelihood_x, loglikelihood_x_i, loglikelihood_x_lessi;
  double logprior_conditional_x, logproposal_conditional_x, log_artificial_likelihood_x;
  double logtarget_x, target_x, loglikelihood_target_x;
  double loglikelihood_gridpoint, loglikelihood_gridpoint_i;
  double logprior_conditional_gridpoint, logproposal_conditional_gridpoint, log_artificial_likelihood_gridpoint;
  double logtarget_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double conditional_mean, conditional_constants, conditional_precision, centered_mean;
  double vi_conditional_constants, vi_conditional_precision, vi_centered_mean;
  
  arma::rowvec xpoint(dimension);
  arma::vec loglikelihood_x_all(dimension);
  arma::vec xpoint_i(1);
  arma::vec gridpoint_i(1);
  arma::vec xpoint_lessi(dimension-1);
  arma::vec dimension_index = arma::regspace(0, dimension-1);
  arma::uvec other_index(dimension-1);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  arma::vec conditional_centered_mean, vi_conditional_centered_mean;
  
  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);

    // Evaluate loglikelihood
    loglikelihood_x = 0;
    for (int i = 0; i < dimension; i++){
      xpoint_i(0) = xpoint(i);
      loglikelihood_x_i = as_scalar(coxprocess_loglikelihood_term(i, xpoint_i, counts));
      loglikelihood_x_all(i) = loglikelihood_x_i;
      loglikelihood_x += loglikelihood_x_i;
    }
    
    for (int i = 0; i < dimension; i++){
      // Set up variables
      xpoint_i(0) = xpoint(i);
      other_index = arma::find(dimension_index != i);
      xpoint_lessi = xpoint.elem(other_index);
      
      // Precompute prior conditional distribution 
      conditional_centered_mean = xpoint_lessi - parameter_mu;
      conditional_mean = parameter_mu + as_scalar(conditional_term.row(i) * conditional_centered_mean);
      conditional_constants = - 0.5 * ( log2pi + log(conditional_var(i)) );
      conditional_precision = 1.0 / conditional_var(i);
      
      // Precompute VI conditional distribution 
      vi_conditional_constants = - 0.5 * ( log2pi + log(vi_var(i)) );
      vi_conditional_precision = 1.0 / vi_var(i);
      
      // Evaluate loglikelihood
      loglikelihood_x_i = loglikelihood_x_all(i);
      loglikelihood_x_lessi = loglikelihood_x - loglikelihood_x_i;
      
      // Evaluate conditional densities
      centered_mean = xpoint(i) - conditional_mean;
      logprior_conditional_x = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
      vi_centered_mean = xpoint(i) - vi_mean(i);
      logproposal_conditional_x = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
      log_artificial_likelihood_x = logprior_conditional_x + loglikelihood_x - logproposal_conditional_x;
      logtarget_x = (1.0 - lambda) * logproposal_conditional_x + lambda * (logprior_conditional_x + loglikelihood_x_i);
      target_x = std::exp(logtarget_x);
      loglikelihood_target_x = log_artificial_likelihood_x * target_x;
      
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
      
      // Evaluate integrands at lower gridpoints
      target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
      target_lowergrid(xindex) = target_x;
      loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
      for (int p = 0; p < xindex; p++){
        gridpoint_i(0) = lowergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = lowergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        vi_centered_mean = lowergrid(p) - vi_mean(i);
        logproposal_conditional_gridpoint = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
      }
      
      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i(0) = uppergrid(p);
        
        // Evaluate loglikelihood
        loglikelihood_gridpoint_i = as_scalar(coxprocess_loglikelihood_term(i, gridpoint_i, counts));
        loglikelihood_gridpoint = loglikelihood_gridpoint_i + loglikelihood_x_lessi;
        
        // Evaluate conditional densities
        centered_mean = uppergrid(p) - conditional_mean;
        logprior_conditional_gridpoint = conditional_constants - 0.5 * conditional_precision * centered_mean * centered_mean;
        vi_centered_mean = uppergrid(p) - vi_mean(i);
        logproposal_conditional_gridpoint = vi_conditional_constants - 0.5 * vi_conditional_precision * vi_centered_mean * vi_centered_mean;
        log_artificial_likelihood_gridpoint = logprior_conditional_gridpoint + loglikelihood_gridpoint - logproposal_conditional_gridpoint;
        logtarget_gridpoint = (1.0 - lambda) * logproposal_conditional_gridpoint + lambda * (logprior_conditional_gridpoint + loglikelihood_gridpoint_i);
        target_gridpoint = std::exp(logtarget_gridpoint);
        target_uppergrid(p) = target_gridpoint;
        loglikelihood_target_uppergrid(p) = log_artificial_likelihood_gridpoint * target_gridpoint;
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
      gibbs_velocity(n, i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
      
    }
  }
  
  return(gibbs_velocity);
  
}

