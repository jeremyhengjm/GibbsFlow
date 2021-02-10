#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Problem settings
const int dimension = 4;
double mixturemodel_double_dimension = static_cast<double>(dimension);
const double mixturemodel_sigmasq = 0.55 * 0.55;
const double mixturemodel_boundary = 10.0;

// Prior
const double mixturemodel_priorconst = - mixturemodel_double_dimension * std::log(20.0);

//' @rdname mixturemodel_logprior
//' @title Evaluate mixture model prior density
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mixturemodel_logprior(arma::mat x){
  int n = x.n_rows;
  arma::vec output(n);

  for (int i = 0; i < n; i++){
    if (all(abs(x.row(i)) < mixturemodel_boundary)){
      output(i) = mixturemodel_priorconst;
    } else {
      output(i) = arma::datum::log_min;
    }
  }

  return(output);
}

//' @rdname mixturemodel_gradlogprior
//' @title Evaluate gradient of mixture model prior density
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat mixturemodel_gradlogprior(arma::mat x){
  arma::mat output;
  output.zeros(size(x));

  return(output);
}

//' @rdname mixturemodel_partialderivative_logprior
//' @title Evaluate partial derivative of mixture model prior density
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double mixturemodel_partialderivative_logprior(int dim, arma::rowvec x){
  return(0);

}

//' @rdname mixturemodel_sampleprior
//' @title Sample from mixture model prior distribution
//' @param n number of samples
//' @return samples
//' @export
// [[Rcpp::export]]
arma::mat mixturemodel_sampleprior(int n){
  arma::mat output(n, dimension);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < dimension; j++){
      output(i, j) = 2.0 * mixturemodel_boundary * arma::randu() - mixturemodel_boundary;
    }
  }

  return(output);
}

// Univariate Gaussian density (in log-scale)
const double normconst = - 0.5 * std::log(2.0 * arma::datum::pi * mixturemodel_sigmasq);
const double normfactor = - 0.5 / mixturemodel_sigmasq;
double lognormpdf(double y, double x){
  double output, yminusx;
  yminusx = y - x;
  output = normconst + normfactor * yminusx * yminusx;

  return(output);
}


// Likelihood

//' @rdname mixturemodel_loglikelihood
//' @title Evaluate mixture model loglikelihood function
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec mixturemodel_loglikelihood(arma::mat x, arma::vec observations, arma::vec mixturelogweights){
  int N = x.n_rows;
  int m = observations.n_elem;
  arma::vec output(N);
  double loglikelihood, maxlogdensities;
  arma::vec logdensities(dimension), densities(dimension);

  for (int n = 0; n < N; n++){
    loglikelihood = 0;
    for (int j = 0; j < m; j++){
      for (int i = 0; i < dimension; i++){
        logdensities(i) = mixturelogweights(i) + lognormpdf(observations(j), x(n, i));
      }
      maxlogdensities = arma::max(logdensities);
      densities = arma::exp(logdensities - maxlogdensities);
      loglikelihood += std::log(sum(densities)) + maxlogdensities;
    }
    output(n) = loglikelihood;
  }

  return(output);
}

//' @rdname mixturemodel_gradloglikelihood
//' @title Evaluate gradient of mixture model loglikelihood function
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat mixturemodel_gradloglikelihood(arma::mat x, arma::vec observations, arma::vec mixturelogweights){
  int N = x.n_rows;
  int m = observations.n_elem;
  arma::mat output = arma::zeros<arma::mat>(N, dimension);
  arma::vec densities(dimension);
  double mixture_sum, current_density;

  for (int n = 0; n < N; n++){
    for (int j = 0; j < m; j++){
      mixture_sum = 0;

      for (int i = 0; i < dimension; i++){
        current_density = std::exp(mixturelogweights(i) + lognormpdf(observations(j), x(n, i)));
        densities(i) = current_density;
        mixture_sum += current_density;
      }

      for (int k = 0; k < dimension; k++){
        output(n, k) += (observations(j) - x(n, k)) * densities(k) / (mixturemodel_sigmasq * mixture_sum);
      }

    }
  }

  return(output);
}

//' @rdname mixturemodel_partialderivative_loglikelihood
//' @title Evaluate partial derivative of mixture model likelihood function
//' @param dim dimension along which partial derivative is taken (dim = 0 is the first coordinate)
//' @param x evaluation point
//' @return partial derivative value
//' @export
// [[Rcpp::export]]
double mixturemodel_partialderivative_loglikelihood(int dim, arma::rowvec x, arma::vec observations, arma::vec mixturelogweights){
  int m = observations.n_elem;
  double output, mixture_sum, current_density, density;
  output = 0;
  for (int j = 0; j < m; j++){
    mixture_sum = 0;
    for (int i = 0; i < dimension; i++){
      current_density = std::exp(mixturelogweights(i) + lognormpdf(observations(j), x(i)));
      mixture_sum += current_density;
      if (i == dim){
        density = current_density;
      }
    }
    output += (observations(j) - x(dim)) * density / (mixturemodel_sigmasq * mixture_sum);
  }

  return(output);
}

// Gibbs flow computation
const double lowerbound = -10.0;
const double upperbound = 10.0;
const int ngridpoints = 100;
const arma::vec grid = arma::linspace(lowerbound, upperbound, ngridpoints);
const double mesh = grid(1) - grid(0);

//' @rdname mixturemodel_gibbsflow
//' @title Compute Gibbs flow for mixture model
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @param observations dataset
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List mixturemodel_gibbsflow(double stepsize, double lambda, double derivative_lambda,
                                arma::mat xparticles, arma::vec logdensity, arma::vec observations, arma::vec mixturelogweights){
  int N = xparticles.n_rows; // no. of particles
  int m = observations.n_elem; // no. of observations

  // Output
  arma::rowvec log_jacobian_dets(N);

  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double log_jacobian_det, maxlogdensities, loglikelihood_x;
  double target_x, loglikelihood_target_x, gridpoint_i, loglikelihood_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;
  double current_density, mixture_sum, density, partialderivative_loglikelihood;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative, new_component;

  arma::rowvec xpoint(dimension);
  arma::mat logdensities(m, dimension);
  arma::rowvec logdensities_row(dimension), densities_row(dimension);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;

  for (int n = 0; n < N; n++){
    xpoint = xparticles.row(n);
    log_jacobian_det = 0;

    // Pre-compute mixture components terms in loglikelihood
    for (int j = 0; j < m; j++){
      for (int i = 0; i < dimension; i++){
        logdensities(j, i) = mixturelogweights(i) + lognormpdf(observations(j), xpoint(i));
      }
    }

    // Compute Gibbs velocity field and move particle
    for (int i = 0; i < dimension; i++){
      // Evaluate log-likelihood and target at particle location
      loglikelihood_x = 0;
      for (int j = 0; j < m; j++){
        logdensities_row = logdensities.row(j);
        maxlogdensities = arma::max(logdensities_row);
        densities_row = arma::exp(logdensities_row - maxlogdensities);
        loglikelihood_x += std::log(sum(densities_row)) + maxlogdensities;
      }
      target_x = std::exp( mixturemodel_priorconst + lambda * loglikelihood_x );
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
        gridpoint_i = lowergrid(p);
        // Evaluate log-likelihood and target at quadrature grid point
        loglikelihood_gridpoint = 0;
        for (int j = 0; j < m; j++){
          logdensities_row = logdensities.row(j);
          logdensities_row(i) = mixturelogweights(i) + lognormpdf(observations(j), gridpoint_i);
          maxlogdensities = arma::max(logdensities_row);
          densities_row = arma::exp(logdensities_row - maxlogdensities);
          loglikelihood_gridpoint += std::log(sum(densities_row)) + maxlogdensities;
        }
        target_gridpoint = std::exp( mixturemodel_priorconst + lambda * loglikelihood_gridpoint );
        target_lowergrid(p) = target_gridpoint;
        loglikelihood_target_lowergrid(p) = loglikelihood_gridpoint * target_gridpoint;
      }

      // Evaluate integrands at upper gridpoints
      target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);

      target_uppergrid(0) = target_x;
      loglikelihood_target_uppergrid(0) = loglikelihood_target_x;

      for (int p = 1; p < nuppergrid; p++){
        gridpoint_i = uppergrid(p);
        // Evaluate log-likelihood and target at quadrature grid point
        loglikelihood_gridpoint = 0;
        for (int j = 0; j < m; j++){
          logdensities_row = logdensities.row(j);
          logdensities_row(i) = mixturelogweights(i) + lognormpdf(observations(j), gridpoint_i);
          maxlogdensities = arma::max(logdensities_row);
          densities_row = arma::exp(logdensities_row - maxlogdensities);
          loglikelihood_gridpoint += std::log(sum(densities_row)) + maxlogdensities;
        }
        target_gridpoint = std::exp( mixturemodel_priorconst + lambda * loglikelihood_gridpoint );
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

      // Compute partial derivative of log-likelihood
      partialderivative_loglikelihood = 0;
      for (int j = 0; j < m; j++){
        mixture_sum = 0;
        for (int k = 0; k < dimension; k++){
          current_density = std::exp(mixturelogweights(k) + logdensities(j, k));
          if (k == i){
            density = current_density;
          }
          mixture_sum += current_density;
        }
        partialderivative_loglikelihood += (observations(j) - xpoint(i)) * density / (mixturemodel_sigmasq * mixture_sum);
      }

      // Compute partial derivative of Gibbs velocity field
      logtarget_partial_derivative = lambda * partialderivative_loglikelihood;
      gibbs_partial_derivative = derivative_lambda * (integral_loglikelihood / integral_target - loglikelihood_x) -
        gibbs_velocity * logtarget_partial_derivative;

      // Jacobian of transport map
      log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));

      // Move component of particle
      new_component = xpoint(i) + stepsize * gibbs_velocity;
      xpoint(i) = new_component;

      // Update mixture components terms in loglikelihood
      for (int j = 0; j < m; j++){
        // logdensities(j, i) = lognormpdf(observations(j), new_component);
        logdensities(j, i) = mixturelogweights(i) + lognormpdf(observations(j), new_component);
      }

    }
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;

  }
  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles,
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);

}

//' @rdname mixturemodel_gibbsvelocity
//' @title Compute Gibbs flow for mixture model
//' @param time time
//' @param xpoint particle position
//' @param exponent exponent of tempering schedule
//' @param observations dataset
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::rowvec mixturemodel_gibbsvelocity(double time, arma::rowvec xpoint, double exponent, arma::vec observations, arma::vec mixturelogweights){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * pow(time, exponent - 1.0);
  
  // Output 
  arma::rowvec gibbs_velocity(dimension);
  int m = observations.n_elem; // no. of observations
  
  // Declare variables
  int xindex, nlowergrid, nuppergrid;
  double maxlogdensities, loglikelihood_x;
  double target_x, loglikelihood_target_x, gridpoint_i, loglikelihood_gridpoint, target_gridpoint;
  double lowerintegral_target, upperintegral_target, integral_target, conditional_cdf;
  double lowerintegral_loglikelihood, upperintegral_loglikelihood, integral_loglikelihood;

  arma::mat logdensities(m, dimension);
  arma::rowvec logdensities_row(dimension), densities_row(dimension);
  arma::vec lowergrid, uppergrid, target_lowergrid, target_uppergrid;
  arma::vec loglikelihood_target_lowergrid, loglikelihood_target_uppergrid;
    
  // Pre-compute mixture components terms in loglikelihood
  for (int j = 0; j < m; j++){
    for (int i = 0; i < dimension; i++){
      logdensities(j, i) = mixturelogweights(i) + lognormpdf(observations(j), xpoint(i));
    }
  }
    
  // Evaluate log-likelihood and target at particle location
  loglikelihood_x = 0;
  for (int j = 0; j < m; j++){
    logdensities_row = logdensities.row(j);
    maxlogdensities = arma::max(logdensities_row);
    densities_row = arma::exp(logdensities_row - maxlogdensities);
    loglikelihood_x += std::log(sum(densities_row)) + maxlogdensities;
  }
  target_x = std::exp( mixturemodel_priorconst + lambda * loglikelihood_x );
  loglikelihood_target_x = loglikelihood_x * target_x;
  
  // Compute Gibbs velocity field and move particle
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
      
    // Evaluate integrands at lower gridpoints
    target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
    loglikelihood_target_lowergrid = arma::zeros<arma::vec>(nlowergrid);
      
    target_lowergrid(xindex) = target_x;
    loglikelihood_target_lowergrid(xindex) = loglikelihood_target_x;
      
    for (int p = 0; p < xindex; p++){
      gridpoint_i = lowergrid(p);
      // Evaluate log-likelihood and target at quadrature grid point
      loglikelihood_gridpoint = 0;
      for (int j = 0; j < m; j++){
        logdensities_row = logdensities.row(j);
        logdensities_row(i) = mixturelogweights(i) + lognormpdf(observations(j), gridpoint_i);
        maxlogdensities = arma::max(logdensities_row);
        densities_row = arma::exp(logdensities_row - maxlogdensities);
        loglikelihood_gridpoint += std::log(sum(densities_row)) + maxlogdensities;
      }
      target_gridpoint = std::exp( mixturemodel_priorconst + lambda * loglikelihood_gridpoint );
      target_lowergrid(p) = target_gridpoint;
      loglikelihood_target_lowergrid(p) = loglikelihood_gridpoint * target_gridpoint;
    }
      
    // Evaluate integrands at upper gridpoints
    target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
    loglikelihood_target_uppergrid = arma::zeros<arma::vec>(nuppergrid);
      
    target_uppergrid(0) = target_x;
    loglikelihood_target_uppergrid(0) = loglikelihood_target_x;
      
    for (int p = 1; p < nuppergrid; p++){
      gridpoint_i = uppergrid(p);
      // Evaluate log-likelihood and target at quadrature grid point
      loglikelihood_gridpoint = 0;
      for (int j = 0; j < m; j++){
        logdensities_row = logdensities.row(j);
        logdensities_row(i) = mixturelogweights(i) + lognormpdf(observations(j), gridpoint_i);
        maxlogdensities = arma::max(logdensities_row);
        densities_row = arma::exp(logdensities_row - maxlogdensities);
        loglikelihood_gridpoint += std::log(sum(densities_row)) + maxlogdensities;
      }
      target_gridpoint = std::exp( mixturemodel_priorconst + lambda * loglikelihood_gridpoint );
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
      
  }
    
  return(gibbs_velocity);
  
}



