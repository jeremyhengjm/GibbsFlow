#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Variance components model on simulated dataset
// Problem settings
const double alpha_0 = -1.0;
const double beta_0 = 2.0;
const double alpha_1 = 4.0; 
const double beta_1 = 4.0; 
const double inv_beta_1 = 1.0 / beta_1;
const double mu_1 = 0.0;
const double sigma_1_sq = 0.01; 
const double sigma_1 = std::sqrt(sigma_1_sq);
const double mu_2 = 0.0;
const double sigma_2_sq = 0.01;
const double sigma_2 = std::sqrt(sigma_2_sq);
const double inv_sigma_1_sq = 1.0 / sigma_1_sq;
const double inv_sigma_2_sq = 1.0 / sigma_2_sq;
const double sigma_e_sq = 0.00434;
const double inv_sigma_e_sq = 1.0 / sigma_e_sq;

// parameters (sigma_{theta}^2, mu, theta) are stored as a vector x
// using the ordering x[0] = sigma_{theta}^2, x[1] = mu, x[2:(dimension-1)] = theta

// Inverse Gamma distribution
double varcomp_invgamma_logdensity(double x, double alpha, double beta){
  double output = alpha * std::log(beta) - std::lgamma(alpha) - (alpha + 1.0) * std::log(x) - beta / x;
  return(output);
}

// Improper inverse Gamma distribution
double varcomp_invgamma_improper_logdensity(double x, double beta){
  double output = - beta / x;
  return(output);
}

// Univariate Gaussian distribution
const double log2pi = std::log(2.0 * arma::datum::pi);
double varcomp_gaussian_logdensity(double x, double mean, double variance){
  double xc, output;
  xc = x - mean;
  output = - 0.5 * xc * xc / variance;
  output += - 0.5 * log2pi - 0.5 * std::log(variance);
  return(output);
}

// Artificial prior
//' @rdname varcomp_artificial_logprior
//' @title Evaluate variance components model artificial prior density on simulated dataset
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec varcomp_artificial_logprior(arma::mat x, int dimension){
  int N = x.n_rows;
  double logdensity, sigma_theta_sq, mu, theta_i;
  arma::vec output(N);

  for (int n = 0; n < N; n++){
    sigma_theta_sq = x(n, 0);
    mu = x(n, 1);
    logdensity = 0;
    logdensity = varcomp_invgamma_logdensity(sigma_theta_sq, alpha_1, beta_1);
    logdensity += varcomp_gaussian_logdensity(mu, mu_1, sigma_1_sq);
    for (int i = 2; i < dimension; i++){
      theta_i = x(n, i);
      logdensity += varcomp_gaussian_logdensity(theta_i, mu_2, sigma_2_sq);
    }
    output(n) = logdensity;

  }
  return(output);
}

//' @rdname varcomp_gradlogprior_artificial
//' @title Evaluate gradient of variance components artificial prior density on simulated dataset
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat varcomp_gradlogprior_artificial(arma::mat x, int dimension){
  int N = x.n_rows;
  double inv_sigma_theta_sq, mu;
  arma::mat output(N, dimension);

  for (int n = 0; n < N; n++){
    inv_sigma_theta_sq = 1.0 / x(n, 0);
    output(n, 0) = (- alpha_1 - 1.0) * inv_sigma_theta_sq + beta_1 * inv_sigma_theta_sq * inv_sigma_theta_sq;
    mu = x(n, 1);
    output(n, 1) = (mu_1 - mu) * inv_sigma_1_sq;
    for (int i = 2; i < dimension; i++){
      output(n, i) = (mu_2 - x(n, i)) * inv_sigma_2_sq;
    }

  }
  return(output);
}

// Prior
//' @rdname varcomp_logprior
//' @title Evaluate variance components model prior density on simulated dataset
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec varcomp_logprior(arma::mat x, int dimension){
  int N = x.n_rows;
  double logdensity, sigma_theta_sq, mu, theta_i;
  arma::vec output(N);

  for (int n = 0; n < N; n++){
    sigma_theta_sq = x(n, 0);
    mu = x(n, 1);
    logdensity = 0;
    logdensity = varcomp_invgamma_improper_logdensity(sigma_theta_sq, beta_0);
    for (int i = 2; i < dimension; i++){
      theta_i = x(n, i);
      logdensity += varcomp_gaussian_logdensity(theta_i, mu, sigma_theta_sq);
    }
    output(n) = logdensity;

  }
  return(output);
}

//' @rdname varcomp_gradlogprior
//' @title Evaluate gradient of variance components prior density on simulated dataset
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat varcomp_gradlogprior(arma::mat x, int dimension){
  int N = x.n_rows;
  int nobservations = dimension - 2;
  double inv_sigma_theta_sq, mu, theta_i_centered;
  double sum_squares, sum_linear;
  arma::mat output(N, dimension);

  for (int n = 0; n < N; n++){
    inv_sigma_theta_sq = 1.0 / x(n, 0);
    mu = x(n, 1);
    sum_squares = 0;
    sum_linear = 0;
    for (int i = 2; i < dimension; i++){
      theta_i_centered = x(n, i) - mu;
      sum_squares += theta_i_centered * theta_i_centered;
      sum_linear += theta_i_centered;
      output(n, i) = - theta_i_centered * inv_sigma_theta_sq;
    }
    output(n, 0) = (- alpha_0 - 1.0 - 0.5 * nobservations) * inv_sigma_theta_sq + (0.5 * sum_squares + beta_0) * inv_sigma_theta_sq * inv_sigma_theta_sq;
    output(n, 1) = sum_linear * inv_sigma_theta_sq;
  }
  return(output);
}

// Sample from artificial prior distribution
//' @rdname varcomp_sample_artificial_prior
//' @title Sample from variance components model artificial prior distribution
//' @param n number of samples
//' @return samples
//' @export
// [[Rcpp::export]]
arma::mat varcomp_sample_artificial_prior(int N, int dimension){
  arma::mat output(N, dimension);
  double sigma_theta_sq, mu;

  for (int n = 0; n < N; n++){
    sigma_theta_sq = 1.0 / R::rgamma(alpha_1, inv_beta_1);
    output(n, 0) = sigma_theta_sq;
    mu = mu_1 + sigma_1 * arma::randn();
    output(n, 1) = mu;
    for (int i = 2; i < dimension; i++){
      output(n, i) = mu_2 + sigma_2 * arma::randn();
    }
  }

  return(output);
}

// Likelihood
//' @rdname varcomp_loglikelihood
//' @title Evaluate variance components model loglikelihood function on simulated dataset
//' @param x evaluation points
//' @return density values
//' @export
// [[Rcpp::export]]
arma::vec varcomp_loglikelihood(arma::mat x, arma::mat dataset){
  int N = x.n_rows;
  int dimension = dataset.n_rows + 2;
  int nrepetitions = dataset.n_cols;

  double loglikelihood, obs_minus_theta;
  arma::vec output(N);
  for (int n = 0; n < N; n++){
    loglikelihood = 0;
    for (int i = 2; i < dimension; i++){
      for (int j = 0; j < nrepetitions; j++){
        obs_minus_theta = dataset(i-2, j) - x(n, i);
        loglikelihood += obs_minus_theta * obs_minus_theta;
      }
    }
    output(n) = - 0.5 * loglikelihood * inv_sigma_e_sq;
  }

  return(output);
}

//' @rdname varcomp_gradloglikelihood
//' @title Evaluate gradient of variance components loglikelihood on simulated dataset
//' @param x evaluation points
//' @return gradient values
//' @export
// [[Rcpp::export]]
arma::mat varcomp_gradloglikelihood(arma::mat x, arma::vec dataset_rowsums, int dimension, int nrepetitions){
  int N = x.n_rows;

  // double theta_i_centered;
  arma::mat output = arma::zeros<arma::mat>(N, dimension);

  for (int n = 0; n < N; n++){
    for (int i = 2; i < dimension; i++){
      output(n, i) = (dataset_rowsums(i-2) - nrepetitions * x(n, i)) * inv_sigma_e_sq;
    }
  }
  return(output);
}

// Gibbs flow computation 
const double lowerbound = arma::datum::eps;
const int ngridpoints = 50;

//' @rdname varcomp_gibbsflow
//' @title Compute Gibbs flow for variance components model on simulated dataset
//' @param stepsize numerical integration step size
//' @param lambda tempering level
//' @param lambda_next next tempering level
//' @param derivative_lambda time derivative of tempering schedule
//' @param xparticles particle locations
//' @param logdensity corresponding log density values
//' @return list with keys:
//' \code{xparticles} new particle locations,
//' \code{log_jacobian_dets} log Jacobian determinant of the flow
//' @export
// [[Rcpp::export]]
Rcpp::List varcomp_gibbsflow(double stepsize, double lambda, double lambda_next, double derivative_lambda,
                              arma::mat xparticles, arma::vec logdensity, arma::vec dataset_rowsums,
                              int dimension){
  int N = xparticles.n_rows; // no. of particles
  int nobservations = dimension - 2;

  // Output
  arma::rowvec log_jacobian_dets(N);

  // Declare variables
  arma::rowvec xpoint;
  double log_jacobian_det;
  double sigma_theta_sq, mu, theta_i, theta_i_centered, sum_squares, sum_thetas;
  double current_alpha, current_beta, d_alpha, d_beta, current_xi, log_sigma_theta_sq;
  double target_x, d_logtarget_x, d_logtarget_target_x;
  double gridpoint, log_gridpoint, target_gridpoint, d_logtarget_gridpoint;
  double gibbs_velocity, logtarget_partial_derivative, gibbs_partial_derivative;
  double combined_var, current_sd, current_mean, next_sd, next_mean, ratio_sd, temp_s_sq, temp_m;

  arma::vec lowergrid(ngridpoints+1);
  arma::vec d_logtarget_target_lowergrid(ngridpoints+1);
  
  // Precompute
  current_alpha = alpha_1 + lambda * (alpha_0 - alpha_1 + 0.5 * nobservations);
  d_alpha = derivative_lambda * (alpha_0 - alpha_1 + 0.5 * nobservations);
  
  for (int n = 0; n < N; n++){
    // Initialize each particle
    xpoint = xparticles.row(n);
    log_jacobian_det = 0;

    // Extract parameters and summary statistics
    sigma_theta_sq = xpoint(0);
    mu = xpoint(1);
    sum_squares = 0;
    sum_thetas = 0;
    for (int i = 2; i < dimension; i++){
      theta_i = xpoint(i);
      theta_i_centered = theta_i - mu;
      sum_squares += theta_i_centered * theta_i_centered;
      sum_thetas += theta_i;
    }

    // Update sigma_{theta}^2 using inverse Gamma flow
    current_beta = beta_1 + lambda * (beta_0 - beta_1 + 0.5 * sum_squares);
    d_beta = derivative_lambda * (beta_0 - beta_1 + 0.5 * sum_squares);
    current_xi = d_alpha * ( std::log(current_beta) - R::digamma(current_alpha) ) +
      d_beta * current_alpha / current_beta;
    log_sigma_theta_sq = std::log(sigma_theta_sq);
    target_x = std::exp( -(current_alpha + 1.0) * log_sigma_theta_sq - current_beta / sigma_theta_sq );
    d_logtarget_x = current_xi - d_alpha * log_sigma_theta_sq - d_beta / sigma_theta_sq;
    d_logtarget_target_x = d_logtarget_x * target_x;

    lowergrid = arma::linspace(lowerbound, sigma_theta_sq, ngridpoints+1);

    // Evaluate integrands at gridpoints on sigma_{theta}^2 space
    d_logtarget_target_lowergrid = arma::zeros<arma::vec>(ngridpoints+1);
    d_logtarget_target_lowergrid(ngridpoints) = d_logtarget_target_x;
    for (int p = 0; p < ngridpoints; p++){
      gridpoint = lowergrid(p);
      log_gridpoint = std::log(gridpoint);
      target_gridpoint = std::exp( -(current_alpha + 1.0) * log_gridpoint - current_beta / gridpoint );
      d_logtarget_gridpoint = current_xi - d_alpha * log_gridpoint - d_beta / gridpoint;
      d_logtarget_target_lowergrid(p) = d_logtarget_gridpoint * target_gridpoint;
    }

    // Compute Gibbs velocity on sigma_{theta}^2 space and move
    gibbs_velocity = - as_scalar(trapz(lowergrid, d_logtarget_target_lowergrid)) / target_x;
    logtarget_partial_derivative = - (current_alpha + 1.0) / sigma_theta_sq + current_beta / (sigma_theta_sq * sigma_theta_sq);
    gibbs_partial_derivative = - d_logtarget_x - gibbs_velocity * logtarget_partial_derivative;
    log_jacobian_det += std::log(std::abs(1.0 + stepsize * gibbs_partial_derivative));
    sigma_theta_sq = sigma_theta_sq + stepsize * gibbs_velocity;
    xpoint(0) = sigma_theta_sq;

    // Update mu using Gaussian flow 
    combined_var = sigma_theta_sq * (1 - lambda) + sigma_1_sq * nobservations * lambda;
    current_sd = std::sqrt(sigma_1_sq * sigma_theta_sq / combined_var);
    current_mean = (sigma_theta_sq * (1 - lambda) * mu_1 + lambda * sigma_1_sq * sum_thetas) / combined_var;

    combined_var = sigma_theta_sq * (1 - lambda_next) + sigma_1_sq * nobservations * lambda_next;
    next_sd = std::sqrt(sigma_1_sq * sigma_theta_sq / combined_var);
    next_mean = (sigma_theta_sq * (1 - lambda_next) * mu_1 + lambda_next * sigma_1_sq * sum_thetas) / combined_var;

    ratio_sd = next_sd / current_sd;
    log_jacobian_det += std::log(ratio_sd);
    mu = ratio_sd * (mu - current_mean) + next_mean;
    xpoint(1) = mu;

    // Update theta's using Gaussian flow
    for (int i = 2; i < dimension; i++){
      theta_i = xpoint(i);
      // Current
      combined_var = lambda * sigma_2_sq + (1 - lambda) * sigma_theta_sq;
      temp_s_sq = sigma_theta_sq * sigma_2_sq / combined_var;
      temp_m = (lambda * mu * sigma_2_sq + (1 - lambda) * mu_2 * sigma_theta_sq) / combined_var;

      combined_var = sigma_e_sq + lambda * temp_s_sq;
      current_sd = std::sqrt(temp_s_sq * sigma_e_sq / combined_var);
      current_mean = (sigma_e_sq * temp_m + lambda * temp_s_sq * dataset_rowsums(i-2)) / combined_var;

      // Next
      combined_var = lambda_next * sigma_2_sq + (1 - lambda_next) * sigma_theta_sq;
      temp_s_sq = sigma_theta_sq * sigma_2_sq / combined_var;
      temp_m = (lambda_next * mu * sigma_2_sq + (1 - lambda_next) * mu_2 * sigma_theta_sq) / combined_var;

      combined_var = sigma_e_sq + lambda_next * temp_s_sq;
      next_sd = std::sqrt(temp_s_sq * sigma_e_sq / combined_var);
      next_mean = (sigma_e_sq * temp_m + lambda_next * temp_s_sq * dataset_rowsums(i-2)) / combined_var;

      ratio_sd = next_sd / current_sd;
      log_jacobian_det += std::log(ratio_sd);
      theta_i = ratio_sd * (theta_i - current_mean) + next_mean;
      xpoint(i) = theta_i;
    }

    // Store new particle position and jacobian determinant
    xparticles.row(n) = xpoint;
    log_jacobian_dets(n) = log_jacobian_det;

  }

  return Rcpp::List::create(Rcpp::Named("xparticles") = xparticles,
                            Rcpp::Named("log_jacobian_dets") = log_jacobian_dets);
}

//' @rdname varcomp_gibbsvelocity
//' @title Compute Gibbs velocity for variance components model on simulated dataset
//' @param time time
//' @param xparticles particle positions
//' @param exponent exponent of tempering schedule
//' @return gibbs_velocity gibbs velocity field
//' @export
// [[Rcpp::export]]
arma::mat varcomp_gibbsvelocity(double time, arma::mat xparticles, double exponent, arma::vec dataset_rowsums, int dimension){
  double lambda = std::pow(time, exponent);
  double derivative_lambda = exponent * std::pow(time, exponent - 1.0);
  int N = xparticles.n_rows; // no. of particles
  int nobservations = dimension - 2;
  
  // Output
  arma::mat gibbs_velocity(N, dimension);
  
  // Declare variables
  double sigma_theta_sq, mu, theta_i, theta_i_centered, sum_squares, sum_thetas;
  double current_alpha, current_beta, d_alpha, d_beta, current_xi, log_sigma_theta_sq;
  double target_x, d_logtarget_x, d_logtarget_target_x;
  double gridpoint, log_gridpoint, target_gridpoint, d_logtarget_gridpoint;
  double combined_var, current_sd, d_current_var, d_current_sd;
  double mean_numerator, mean_denominator, current_mean, d_mean_numerator, d_mean_denominator, d_mean;  
  double m_numerator, m_denominator, m_summary, d_m_numerator, d_m_denominator, d_m_summary;
  double s_sq, d_s_sq;
  double tau_sq_numerator, tau_sq_denominator, tau_sq, tau, d_tau_sq_numerator, d_tau_sq_denominator, d_tau_sq, d_tau;
  double xi_numerator, xi_denominator, xi_summary, d_xi_numerator, d_xi_denominator, d_xi_summary;
  
  arma::vec lowergrid(ngridpoints+1);
  arma::vec d_logtarget_target_lowergrid(ngridpoints+1);
  arma::rowvec xpoint;
  
  // precompute
  current_alpha = alpha_1 + lambda * (alpha_0 - alpha_1 + 0.5 * nobservations);
  d_alpha = derivative_lambda * (alpha_0 - alpha_1 + 0.5 * nobservations);
  
  for (int n = 0; n < N; n++){
    // Initialize each particle
    xpoint = xparticles.row(n);
    
    // Extract parameters and summary statistics
    sigma_theta_sq = xpoint(0);
    mu = xpoint(1);
    sum_squares = 0;
    sum_thetas = 0;
    for (int i = 2; i < dimension; i++){
      theta_i = xpoint(i);
      theta_i_centered = theta_i - mu;
      sum_squares += theta_i_centered * theta_i_centered;
      sum_thetas += theta_i;
    }
    
    // Update sigma_{theta}^2 using inverse Gamma flow
    current_beta = beta_1 + lambda * (beta_0 - beta_1 + 0.5 * sum_squares);
    d_beta = derivative_lambda * (beta_0 - beta_1 + 0.5 * sum_squares);
    current_xi = d_alpha * ( std::log(current_beta) - R::digamma(current_alpha) ) +
      d_beta * current_alpha / current_beta;
    log_sigma_theta_sq = std::log(sigma_theta_sq);
    target_x = std::exp( -(current_alpha + 1.0) * log_sigma_theta_sq - current_beta / sigma_theta_sq );
    d_logtarget_x = current_xi - d_alpha * log_sigma_theta_sq - d_beta / sigma_theta_sq;
    d_logtarget_target_x = d_logtarget_x * target_x;
    
    lowergrid = arma::linspace(lowerbound, sigma_theta_sq, ngridpoints+1);
    
    // Evaluate integrands at gridpoints on sigma_{theta}^2 space
    d_logtarget_target_lowergrid = arma::zeros<arma::vec>(ngridpoints+1);
    d_logtarget_target_lowergrid(ngridpoints) = d_logtarget_target_x;
    for (int p = 0; p < ngridpoints; p++){
      gridpoint = lowergrid(p);
      log_gridpoint = std::log(gridpoint);
      target_gridpoint = std::exp( -(current_alpha + 1.0) * log_gridpoint - current_beta / gridpoint );
      d_logtarget_gridpoint = current_xi - d_alpha * log_gridpoint - d_beta / gridpoint;
      d_logtarget_target_lowergrid(p) = d_logtarget_gridpoint * target_gridpoint;
    }
    
    // Compute Gibbs velocity in sigma_{theta}^2 space 
    gibbs_velocity(n, 0) = - as_scalar(trapz(lowergrid, d_logtarget_target_lowergrid)) / target_x;
    
    // Update mu using Gaussian flow
    combined_var = sigma_theta_sq * (1.0 - lambda) + sigma_1_sq * nobservations * lambda;
    current_sd = std::sqrt(sigma_1_sq * sigma_theta_sq / combined_var);
    d_current_var = - sigma_1_sq * sigma_theta_sq * (- sigma_theta_sq * derivative_lambda + sigma_1_sq * nobservations * derivative_lambda) / (combined_var * combined_var);
    d_current_sd = 0.5 * d_current_var / current_sd; // chain rule
    
    mean_numerator = sigma_theta_sq * (1.0 - lambda) * mu_1 + lambda * sigma_1_sq * sum_thetas;
    mean_denominator = combined_var;
    current_mean = mean_numerator / mean_denominator;
    d_mean_numerator = - sigma_theta_sq * derivative_lambda * mu_1 + derivative_lambda * sigma_1_sq * sum_thetas;
    d_mean_denominator = - sigma_theta_sq * derivative_lambda + sigma_1_sq * nobservations * derivative_lambda;
    d_mean = (d_mean_numerator * mean_denominator - mean_numerator * d_mean_denominator) / (mean_denominator * mean_denominator); // quotient rule 
    
    // Compute Gibbs velocity in mu space
    gibbs_velocity(n, 1) = d_current_sd * (mu - current_mean) / current_sd + d_mean;
    
    // Update theta's using Gaussian flow
    for (int i = 2; i < dimension; i++){
      theta_i = xpoint(i);
      
      // m summary 
      m_numerator = lambda * mu * sigma_2_sq + (1.0 - lambda) * mu_2 * sigma_theta_sq;
      m_denominator = lambda * sigma_2_sq + (1.0 - lambda) * sigma_theta_sq;
      m_summary = m_numerator / m_denominator;
      
      d_m_numerator = derivative_lambda * mu * sigma_2_sq - derivative_lambda * mu_2 * sigma_theta_sq;
      d_m_denominator = derivative_lambda * sigma_2_sq - derivative_lambda * sigma_theta_sq;
      d_m_summary = (d_m_numerator * m_denominator - m_numerator * d_m_denominator) / (m_denominator * m_denominator); // quotient rule 
      
      // s_sq summary  
      s_sq = sigma_theta_sq * sigma_2_sq / m_denominator;
      d_s_sq = - sigma_theta_sq * sigma_2_sq * (derivative_lambda * sigma_2_sq - derivative_lambda * sigma_theta_sq) / (m_denominator * m_denominator);
      
      // tau summary
      tau_sq_numerator = s_sq * sigma_e_sq;
      tau_sq_denominator = sigma_e_sq + lambda * s_sq;
      tau_sq = tau_sq_numerator / tau_sq_denominator; 
      tau = std::sqrt(tau_sq); 
      
      d_tau_sq_numerator = d_s_sq * sigma_e_sq;
      d_tau_sq_denominator = derivative_lambda * s_sq + lambda * d_s_sq;
      d_tau_sq = (d_tau_sq_numerator * tau_sq_denominator - tau_sq_numerator * d_tau_sq_denominator) / (tau_sq_denominator * tau_sq_denominator); // quotient rule 
      d_tau = 0.5 * d_tau_sq / tau; // chain rule
      
      // xi summary
      xi_numerator = sigma_e_sq * m_summary + lambda * s_sq * dataset_rowsums(i-2);
      xi_denominator = sigma_e_sq + lambda * s_sq;   
      xi_summary = xi_numerator / xi_denominator;
      
      d_xi_numerator = sigma_e_sq * d_m_summary + (derivative_lambda * s_sq + lambda * d_s_sq) * dataset_rowsums(i-2);
      d_xi_denominator = derivative_lambda * s_sq + lambda * d_s_sq;   
      d_xi_summary = (d_xi_numerator * xi_denominator - xi_numerator * d_xi_denominator) / (xi_denominator * xi_denominator); // quotient rule
      
      // Compute Gibbs velocity in theta_i space
      gibbs_velocity(n, i) = d_tau * (theta_i - xi_summary) / tau + d_xi_summary;
      
    } 
    
  }
  
  return(gibbs_velocity);  
  
}

