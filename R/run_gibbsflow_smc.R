#' @rdname run_gibbsflow_smc
#' @title Run Gibbs flow sequential Monte Carlo sampler
#' @param prior list with keys: 
#' \code{logdensity} evaluates log prior density, 
#' \code{gradlogdensity} returns its gradient, 
#' \code{rinit} samples from the prior distribution
#' @param likelihood list with keys: 
#' \code{logdensity} samples from proposal, 
#' \code{gradlogdensity} returns its gradient
#' @param nparticles number of particles
#' @param timegrid vector describing numerical integration times 
#' @param lambda vector describing tempering schedule
#' @param derivative_lambda time derivative of tempering schedule
#' @param compute_gibbsflow function computing Gibbs flow
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
#' @seealso \code{\link{run_gibbsflow_ais}} if no resampling is desired
#' @export
run_gibbsflow_smc <- function(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, mcmc){
  # initialization
  xparticles <- prior$rinit(nparticles)
  previous_logdensity <- prior$logdensity(xparticles)
  
  # pre-allocate
  dimension <- ncol(xparticles)
  nsteps <- length(lambda) # same length as timegrid
  stepsize <- diff(timegrid) 
  # xtrajectory <- array(dim = c(nparticles, dimension, nsteps))
  # xtrajectory[ , , 1] <- xparticles
  logweights <- rep(0, nparticles)
  ess <- rep(0, nsteps)
  ess[1] <- nparticles
  log_normconst <- rep(0, nsteps)
  log_ratio_normconst <- 0
  acceptprob <- matrix(nrow = 2, ncol = nsteps)
  acceptprob[ , 1] <- c(1, 1)
  
  for (istep in 2:nsteps){
    # gibbs flow move
    output_flow <- compute_gibbsflow(stepsize[istep-1], lambda[istep-1], lambda[istep], derivative_lambda[istep-1], 
                                     xparticles, previous_logdensity)
    xparticles <- output_flow$xparticles 
    log_jacobian_dets <- as.numeric(output_flow$log_jacobian_dets)
    
    # weight 
    current_logdensity <- prior$logdensity(xparticles) + lambda[istep] * likelihood$logdensity(xparticles)
    logweights <- current_logdensity - previous_logdensity + log_jacobian_dets
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
    current_logdensity <- current_logdensity[ancestors]

    # MCMC
    current_logtarget <- function(x) prior$logdensity(x) + lambda[istep] * likelihood$logdensity(x)
    current_gradlogtarget <- function(x) prior$gradlogdensity(x) + lambda[istep] * likelihood$gradlogdensity(x)
    transition_kernel <- construct_kernel(current_logtarget, current_gradlogtarget, mcmc)
    
    current_acceptprob <- rep(0, mcmc$nmoves)
    for (imove in 1:mcmc$nmoves){
      mcmc_output <- transition_kernel(xparticles, current_logdensity)
      xparticles <- mcmc_output$x
      current_logdensity <- mcmc_output$logtarget_x
      current_acceptprob[imove] <- mcmc_output$acceptprob
    }
    acceptprob[ , istep] <- c(min(current_acceptprob), max(current_acceptprob))
    
    # store trajectory
    # xtrajectory[ , , istep] <- xparticles
    previous_logdensity <- current_logdensity
    
  }
  
  return(list(xparticles = xparticles, ess = ess, log_normconst = log_normconst, acceptprob = acceptprob))
}

