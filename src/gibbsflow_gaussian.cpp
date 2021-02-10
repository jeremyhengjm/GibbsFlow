#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Problem settings
const int dimension = 2;
  
// Prior 
const double log2pi = std::log(2.0 * arma::datum::pi);
const arma::rowvec prior_mean = arma::zeros<arma::rowvec>(dimension);
const arma::mat prior_cov = arma::eye(dimension, dimension);
const arma::mat prior_precision = arma::inv_sympd(prior_cov);
const arma::mat prior_chol = arma::chol(prior_cov);
const arma::mat prior_rooti = arma::trans(arma::inv(arma::trimatu(prior_chol)));
const double prior_rootisum = arma::sum(log(prior_rooti.diag()));

//' @rdname gaussian_logprior
//' @title Evaluate multivariate Gaussian prior density
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec gaussian_logprior(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = -(static_cast<double>(dimension)/2.0) * log2pi + prior_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = prior_rooti * arma::trans(x.row(i) - prior_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname gaussian_gradlogprior
//' @title Evaluate gradient of multivariate Gaussian prior density
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat gaussian_gradlogprior(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);

  for (int i = 0; i < n; i++){
    output.row(i) = (prior_mean - x.row(i)) * prior_precision;
  }
  
  return(output);
}

//' @rdname gaussian_partialderivative_logprior
//' @title Evaluate partial derivative of multivariate Gaussian prior density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double gaussian_partialderivative_logprior(int dim, arma::rowvec x){
  arma::rowvec xc = prior_mean - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += prior_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname gaussian_sampleprior
//' @title Sample from multivariate Gaussian prior distribution
//' @param n number of samples 
//' @return samples 
//' @export
// [[Rcpp::export]]
arma::mat gaussian_sampleprior(int n){
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    arma::rowvec Z = arma::randn(1, dimension);
    output.row(i) = prior_mean + Z * prior_chol;
  }
  
  return(output);
}

// Likelihood
// const arma::rowvec like_obs = 0 * arma::ones(1, dimension);
const arma::rowvec like_obs = 14.25 * arma::ones(1, dimension);
const double like_rho = 0.5;
const double like_var = 1.0;
arma::mat construct_like_cov(int dimension){
  arma::mat like_cov = arma::zeros<arma::mat>(dimension, dimension);
  for (int i = 0; i < dimension; i++){
    for (int j = 0; j < dimension; j++){
      if (i == j) {
        like_cov(i,j) = like_var;
      } else {
        like_cov(i,j) = like_rho;
      }
    }
  }
  return(like_cov);
}
const arma::mat like_cov = construct_like_cov(dimension);
const arma::mat like_precision = arma::inv(like_cov);
const arma::mat like_rooti = arma::trans(arma::inv(arma::trimatu(arma::chol(like_cov))));
const double like_rootisum = arma::sum(log(like_rooti.diag()));

//' @rdname gaussian_loglikelihood
//' @title Evaluate multivariate Gaussian loglikelihood function
//' @param x evaluation points  
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec gaussian_loglikelihood(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);
  double constants = -(static_cast<double>(dimension)/2.0) * log2pi + like_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = like_rooti * arma::trans(x.row(i) - like_obs);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
}

//' @rdname gaussian_gradloglikelihood
//' @title Evaluate gradient of multivariate Gaussian loglikelihood function
//' @param x evaluation points  
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat gaussian_gradloglikelihood(arma::mat x){
  int n = x.n_rows;
  arma::mat output(n, dimension);
  
  for (int i = 0; i < n; i++){
    output.row(i) = (like_obs - x.row(i)) * like_precision;
  }
  
  return(output);
}

//' @rdname gaussian_partialderivative_loglikelihood
//' @title Evaluate partial derivative of multivariate Gaussian likelihood function
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point  
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double gaussian_partialderivative_loglikelihood(int dim, arma::rowvec x){
  arma::rowvec xc = like_obs - x;
  double output = 0;
  
  for (int i = 0; i < dimension; i++){
    output += like_precision(dim, i) * xc(i);
  }
  
  return(output);
}

//' @rdname gaussian_posterior_cov
//' @title Compute posterior covariance of multivariate Gaussian example
//' @param lambda temperature
//' @return posterior covariance
//' @export
// [[Rcpp::export]]
arma::mat gaussian_posterior_cov(double lambda){
  arma::mat output = arma::inv_sympd(prior_precision + lambda * like_precision);
  return(output);
}

//' @rdname gaussian_posterior_mean
//' @title Compute posterior mean of multivariate Gaussian example
//' @param lambda temperature
//' @return posterior mean
//' @export
// [[Rcpp::export]]
arma::rowvec gaussian_posterior_mean(double lambda){
  arma::rowvec output = (prior_mean * prior_precision + lambda * like_obs * like_precision) * gaussian_posterior_cov(lambda);
  return(output);
}

//' @rdname gaussian_log_normconst
//' @title Compute log normalizing constant of multivariate Gaussian example
//' @param lambda temperature
//' @return log normalizing constant
//' @export
// [[Rcpp::export]]
double gaussian_log_normconst(double lambda){
  arma::rowvec posterior_mean = gaussian_posterior_mean(lambda);
  arma::mat posterior_cov = gaussian_posterior_cov(lambda);
  arma::mat posterior_precision = prior_precision + lambda * like_precision;
  
  double logdet_priorcov = std::log(arma::det(prior_cov));
  double logdet_likecov = std::log(arma::det(like_cov));
  double logdet_posteriorcov = std::log(arma::det(posterior_cov));
  
  double output = - lambda * ( 0.5 * static_cast<double>(dimension) * log2pi + 0.5 * logdet_likecov) -
    0.5 * logdet_priorcov + 0.5 * logdet_posteriorcov - 
    0.5 * arma::sum( prior_mean % (prior_mean * prior_precision) ) - 
    0.5 * lambda * arma::sum( like_obs %  (like_obs * like_precision) ) + 
    0.5 * sum( posterior_mean % (posterior_mean * posterior_precision) );
  
  return(output);
}

//' @rdname gaussian_logposterior
//' @title Evaluate multivariate Gaussian posterior density
//' @param x evaluation points  
//' @param lambda temperature
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec gaussian_logposterior(arma::mat x, double lambda){
  arma::rowvec posterior_mean = gaussian_posterior_mean(lambda);
  arma::mat posterior_precision = prior_precision + lambda * like_precision;
  arma::mat posterior_cov = arma::inv_sympd(posterior_precision);
  arma::mat posterior_chol = arma::chol(posterior_cov);
  arma::mat posterior_rooti = arma::trans(arma::inv(arma::trimatu(posterior_chol)));
  double posterior_rootisum = arma::sum(log(posterior_rooti.diag()));
  
  int n = x.n_rows;
  arma::vec output(n);
  double constants = -(static_cast<double>(dimension)/2.0) * log2pi + posterior_rootisum;
  
  for (int i = 0; i < n; i++){
    arma::vec z = posterior_rooti * arma::trans(x.row(i) - posterior_mean);
    output(i) = constants - 0.5 * arma::sum(z % z);
  }
  
  return(output);
  
}

// Gibbs flow computation
const double lowerbound = -20.0;
const double upperbound = 20.0;
const int ngridpoints = 500;
const arma::vec grid = arma::linspace(lowerbound, upperbound, ngridpoints);
const double mesh = grid(1) - grid(0);

//' @rdname gaussian_gibbsflow
//' @title Compute Gibbs flow for Gaussian model
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
Rcpp::List gaussian_gibbsflow(double stepsize, double lambda, double derivative_lambda, 
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
      loglikelihood_x = as_scalar(gaussian_loglikelihood(xpoint));
      target_x = std::exp( as_scalar(gaussian_logprior(xpoint)) + lambda * loglikelihood_x );
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
        loglikelihood_gridpoint = as_scalar(gaussian_loglikelihood(gridpoint));
        target_gridpoint = std::exp( as_scalar(gaussian_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
        loglikelihood_gridpoint = as_scalar(gaussian_loglikelihood(gridpoint));
        target_gridpoint = std::exp( as_scalar(gaussian_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
      // gibbs_velocity = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / (target_x + arma::datum::eps); // add machine epsilon
      logtarget_partial_derivative = gaussian_partialderivative_logprior(i, xpoint) + lambda * gaussian_partialderivative_loglikelihood(i, xpoint);
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - loglikelihood_x) - 
                                    gibbs_velocity * logtarget_partial_derivative;
      // gibbs_partial_derivative = gibbs_partial_derivative * target_x / (target_x + arma::datum::eps); // add machine epsilon
      // Rcpp::Rcout << "factor:" << std::endl << target_x / (target_x + arma::datum::eps) << std::endl;
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
      xpoint(i) = xpoint(i) + stepsize * gibbs_velocity;
      
    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;
    
  }
  
  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles, 
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);
  
}

//' @rdname gaussian_gibbsvelocity 
//' @title Compute Gibbs velocity for Gaussian model
//' @param time time
//' @param xpoint particle position
//' @param exponent exponent of tempering schedule
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::rowvec gaussian_gibbsvelocity(double time, arma::rowvec xpoint, double exponent){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * std::pow(time, exponent - 1.0);
  
  // Output 
  arma::rowvec gibbs_velocity(dimension);
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double loglikelihood_x, target_x, loglikelihood_target_x, loglikelihood_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  
  arma::rowvec gridpoint(dimension);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
  
  loglikelihood_x = as_scalar(gaussian_loglikelihood(xpoint));
  target_x = std::exp( as_scalar(gaussian_logprior(xpoint)) + lambda * loglikelihood_x );
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
      loglikelihood_gridpoint = as_scalar(gaussian_loglikelihood(gridpoint));
      target_gridpoint = std::exp( as_scalar(gaussian_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
      loglikelihood_gridpoint = as_scalar(gaussian_loglikelihood(gridpoint));
      target_gridpoint = std::exp( as_scalar(gaussian_logprior(gridpoint)) + lambda * loglikelihood_gridpoint );
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
    
    // Compute Gibbs velocity field 
    gibbs_velocity(i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / target_x;
    // gibbs_velocity(i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / (target_x + std::pow(10, -10));
    // gibbs_velocity(i) = derivative_lambda * (conditional_cdf * integral_loglikelihood - lowerintegral_loglikelihood) / (target_x + arma::datum::eps);
    
  }
  
  return(gibbs_velocity);
  
}

//' @rdname gaussian_stein_variational_importance_sampling
//' @title Perform Stein variational importance sampling
//' @param nfparticles number of follower particles
//' @param nlparticles number of leader particles
//' @param niterations number of time iterations 
//' @param stepsize stepsize of each iteration
//' @return list with keys:
//' \code{xparticles} new particle locations
//' \code{ess} effective sample size
//' \code{log_normconst} log normalizing constant estimate
//' @export
// [[Rcpp::export]]
Rcpp::List gaussian_stein_variational_importance_sampling(int nfparticles, int nlparticles, int niterations, double stepsize, double bandwidth){
  int counter;
  double double_nlparticles = static_cast<double>(nlparticles);
  double difference, distance;//, bandwidth;
  double drift, partial_logtarget, kernel_value, grad_kernel, log_jacobian_det, log_jacobian_det_temp;
  double maxlogweights, ess, log_normconst;
  arma::mat dist_xl_xl(nlparticles, nlparticles);
  arma::vec dist_xl_xl_vec(nlparticles * nlparticles);
  arma::mat dist_xl_xf(nlparticles, nfparticles);
  arma::mat xlparticles = gaussian_sampleprior(nlparticles);
  arma::mat xfparticles = gaussian_sampleprior(nfparticles);
  arma::vec logproposal = gaussian_logprior(xfparticles);
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
    // bandwidth = 0.5 * std::pow(arma::median(dist_xl_xl_vec), 2) / std::log(double_nlparticles + 1.0);
    
    // Move follower particles
    for (int n = 0; n < nfparticles; n++){
      log_jacobian_det = 0;
      
      for (int i = 0; i < dimension; i++){
        drift = 0;
        log_jacobian_det_temp = 0;
        for (int j = 0; j < nlparticles; j++){
          partial_logtarget = gaussian_partialderivative_logprior(i, xlparticles.row(j)) +
            gaussian_partialderivative_loglikelihood(i, xlparticles.row(j));
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
      logweights(n) = as_scalar(gaussian_logprior(xfparticles.row(n))) + as_scalar(gaussian_loglikelihood(xfparticles.row(n))) - logproposal(n);
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
          partial_logtarget = gaussian_partialderivative_logprior(i, xlparticles.row(j)) +
            gaussian_partialderivative_loglikelihood(i, xlparticles.row(j));
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
