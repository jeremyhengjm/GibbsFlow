gaussian_model <- function(dimension){
  
  # construct standard Gaussian prior distribution
  prior_mean <- rep(0, dimension)
  prior_cov <- diag(dimension)
  prior_precision <- solve(prior_cov)
  logprior <- function(x){
    return(mvnpdf(x, prior_mean, prior_cov))  
  }
  gradlogprior <- function(x){
    return(prior_mean - x)
  }
  sampleprior <- function(nparticles){
    return(mvnrnd(nparticles, prior_mean, prior_cov))
  }
  prior <- list(logdensity = logprior, gradlogdensity = gradlogprior, rinit = sampleprior)
  
  # construct Gaussian log-likelihood
  like_obs <- rep(14.25, dimension)
  like_rho <- 0.5
  like_cov <- like_rho * matrix(1, dimension, dimension) + diag(1 - like_rho, dimension)
  like_rooti <- t(solve(chol(like_cov)))
  like_rootisum <- sum(log(diag(like_rooti)))
  loglikelihood <- function(x){
    return(mvnpdf_chol(x, like_obs, like_rooti, like_rootisum))
  }
  like_precision <- solve(like_cov)
  gradloglikelihood <- function(x){
    return((like_obs-x) %*% like_precision)
  }
  likelihood <- list(logdensity = loglikelihood, gradlogdensity = gradloglikelihood)
  
  # construct Gaussian posterior distribution
  posterior <- list()
  posterior_cov <- function(lambda) solve(prior_precision + lambda * like_precision)
  posterior_mean <- function(lambda) (prior_mean %*% prior_precision + lambda * like_obs %*% like_precision) %*% posterior_cov(lambda)
  posterior_log_normconst <- function(lambda){
    current_mean <- posterior_mean(lambda)
    current_precision <- solve(posterior_cov(lambda))
    output <- - lambda * (- 0.5 * dimension * log(2 * pi) + like_rootisum ) -
      0.5 * log(det(prior_cov)) + 0.5 * log(det(posterior_cov(lambda))) - 
      0.5 * sum(prior_mean * (prior_mean %*% prior_precision)) - 
      0.5 * lambda * sum(like_obs *  (like_obs %*% like_precision)) + 
      0.5 * sum(current_mean *  (current_mean %*% current_precision)) 
    return(output)
  }
  posterior_logdensity <- function(x, lambda){
    return(mvnpdf(x, posterior_mean(lambda), posterior_cov(lambda)))  
  } 
  posterior <- list(mean = posterior_mean, cov = posterior_cov, logdensity = posterior_logdensity, 
                    log_normconst = posterior_log_normconst)
  gaussian_model <- list(prior = prior, likelihood = likelihood, posterior = posterior)
  
  return(gaussian_model)
}
