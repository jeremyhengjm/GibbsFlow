#' @rdname run_gibbsflow_sisr
#' @title Run Gibbs flow importance sampler with resampling
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
#' @return list with keys: 
#' \code{xtrajectory} trajectories, 
#' \code{xparticles} particles at terminal time,
#' \code{ess} effective sample size, 
#' \code{log_normconst} log normalizing constant, 
#' @seealso \code{\link{run_gibbsflow_sis}} if no resampling is desired
#' @export
run_gibbsflow_sisr <- function(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow){
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
    previous_logdensity <- current_logdensity[ancestors]
    
    # store trajectory
    # xtrajectory[ , , istep] <- xparticles

  }
  
  return(list(xparticles = xparticles, ess = ess, log_normconst = log_normconst))
}

