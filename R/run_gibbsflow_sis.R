#' @rdname run_gibbsflow_sis 
#' @title Run Gibbs flow sequential importance sampler
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
#' @seealso \code{\link{run_gibbsflow_smc}} if resampling is desired
#' @export
run_gibbsflow_sis <- function(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, gibbsvelocity = NULL){
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
  
  # store norm of velocity
  if (!is.null(gibbsvelocity)){
    normvelocity <- matrix(0, nparticles, nsteps)
    normvelocity[, 1] <- sqrt(rowSums(gibbsvelocity(0, xparticles)^2))
  }
  
  for (istep in 2:nsteps){
    # cat("Time step:", istep, "/", nsteps, "\n")
    # gibbs flow move
    output_flow <- compute_gibbsflow(stepsize[istep-1], lambda[istep-1], lambda[istep], derivative_lambda[istep-1], 
                                     xparticles, previous_logdensity)
    xparticles <- output_flow$xparticles 
    log_jacobian_dets <- as.numeric(output_flow$log_jacobian_dets)

    # weight 
    current_logdensity <- prior$logdensity(xparticles) + lambda[istep] * likelihood$logdensity(xparticles)
    logweights <- logweights + current_logdensity - previous_logdensity + log_jacobian_dets
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)  
    
    # compute effective sample size
    ess[istep] <- 1 / sum(normweights^2)
    
    # compute normalizing constant
    log_normconst[istep] <- log(mean(weights)) + maxlogweights
    
    # store trajectory
    # xtrajectory[ , , istep] <- xparticles
    previous_logdensity <- current_logdensity
    
    # store norm of velocity
    if (!is.null(gibbsvelocity)){
      normvelocity[, istep] <- sqrt(rowSums(gibbsvelocity(timegrid[istep], xparticles)^2))
    }
    
  }
  
  # output
  if (!is.null(gibbsvelocity)){
    return(list(xparticles = xparticles, ess = ess, log_normconst = log_normconst, normvelocity = normvelocity))
  } else {
    return(list(xparticles = xparticles, ess = ess, log_normconst = log_normconst))
  }

}

