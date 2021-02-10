#' @rdname run_smc 
#' @title Run sequential Monte Carlo sampler
#' @param prior list with keys: 
#' \code{logdensity} evaluates log prior density, 
#' \code{gradlogdensity} returns its gradient, 
#' \code{rinit} samples from the prior distribution
#' @param likelihood list with keys: 
#' \code{logdensity} samples from proposal, 
#' \code{gradlogdensity} returns its gradient
#' @param nparticles number of particles
#' @param lambda vector describing tempering schedule
#' @param mcmc list with keys: 
#' \code{choice} specifies type of MCMC method, 
#' \code{parameters} specifies algorithmic tuning parameters,
#' \code{nmoves} specifies number of MCMC move per temperature
#' @return list with keys: 
#' \code{xtrajectory} trajectories, 
#' \code{xparticles} particles at terminal time,
#' \code{ess} effective sample size, 
#' \code{log_normconst} log normalizing constant, 
#' \code{acceptprob} MCMC acceptance probabilities 
#' @seealso \code{\link{run_ais}} if no resampling is desired
#' @export
run_smc <- function(prior, likelihood, nparticles, lambda, mcmc){
  # initialization
  xparticles <- prior$rinit(nparticles)
  
  # pre-allocate
  dimension <- ncol(xparticles)
  nsteps <- length(lambda)
  # xtrajectory <- array(dim = c(nparticles, dimension, nsteps))
  # xtrajectory[ , , 1] <- xparticles
  ess <- rep(0, nsteps)
  ess[1] <- nparticles
  log_normconst <- rep(0, nsteps)
  log_ratio_normconst <- 0
  acceptprob <- matrix(nrow = 2, ncol = nsteps)
  acceptprob[ , 1] <- c(1, 1)
  
  for (istep in 2:nsteps){
    cat("Iteration:", istep, "\n")
    
    # weight 
    current_loglikelihood <- likelihood$logdensity(xparticles)
    logweights <- (lambda[istep] - lambda[istep - 1]) * current_loglikelihood
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    
    # compute effective sample size
    ess[istep] <- 1 / sum(normweights^2)
    
    # compute normalizing constant
    log_ratio_normconst <- log_ratio_normconst + log(mean(weights)) + maxlogweights  
    log_normconst[istep] <- log_ratio_normconst
    
    # resample
    ancestors <- systematic_resampling(normweights, nparticles, runif(1) / nparticles)
    xparticles <- xparticles[ancestors, ]
    current_loglikelihood <- current_loglikelihood[ancestors]

    # MCMC
    current_logtarget <- function(x) prior$logdensity(x) + lambda[istep] * likelihood$logdensity(x)
    current_gradlogtarget <- function(x) prior$gradlogdensity(x) + lambda[istep] * likelihood$gradlogdensity(x)
    transition_kernel <- construct_kernel(current_logtarget, current_gradlogtarget, mcmc)
    
    current_logdensity <- prior$logdensity(xparticles) + lambda[istep] * current_loglikelihood
    current_acceptprob <- rep(0, mcmc$nmoves)
    for (imove in 1:mcmc$nmoves){
      mcmc_output <- transition_kernel(xparticles, current_logdensity)
      xparticles <- mcmc_output$x
      current_logdensity <- mcmc_output$logtarget_x
      current_acceptprob[imove] <- mcmc_output$acceptprob
    }
    acceptprob[ , istep] <- c(min(current_acceptprob), max(current_acceptprob))
    # xtrajectory[ , , istep] <- xparticles
    
  }
  
  return(list(xparticles = xparticles, ess = ess, log_normconst = log_normconst, acceptprob = acceptprob))
}

