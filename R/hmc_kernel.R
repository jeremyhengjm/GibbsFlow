#' @rdname get_hmc_kernel
#' @title Construct Hamiltonian Monte Carlo kernel
#' @param logtarget evaluates log density of the target distribution  
#' @param gradlogtarget evaluates gradient of the log target density 
#' @param parameters list with keys: 
#' \code{stepsize} specifies stepsize of leap-frog discretization, 
#' \code{nsteps} specifies number of leap-frog steps
#' @return list with key: \code{kernel}
#' @export
get_hmc_kernel <- function(logtarget, gradlogtarget, parameters){
  stepsize <- parameters$stepsize
  nsteps <- parameters$nsteps
  
  # leap frog integrator
  leapfrog <- function(x, v){
    # x: N x d matrix; v: N x d matrix
    v <- v + stepsize * gradlogtarget(x) / 2
    for (istep in 1:nsteps){
      x <- x + stepsize * v
      if (istep != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    return(list(x = x, v = v))
  }
  
  # One step of HMC
  kernel <- function(x, logtarget_x){
    # x: N x d matrix
    nparticles <- nrow(x)
    dimension <- ncol(x)
    
    # sample velocity
    v <- matrix(rnorm(nparticles * dimension), nrow = nparticles, ncol = dimension) # could be faster?
    
    # run leap-frog integrator
    leapfrog_result <- leapfrog(x, v)
    proposed_x <- leapfrog_result$x
    proposed_v <- leapfrog_result$v
    
    # compute acceptance probability
    logtarget_proposal <- logtarget(proposed_x)
    accept_ratio <- logtarget_proposal - logtarget_x + rowSums(v^2) / 2 - rowSums(proposed_v^2) / 2
    finite_logical <- is.finite(accept_ratio) # if there is stability issues
    nfinite <- sum(finite_logical) # if there is stability issues
    
    # accept or reject proposal
    accept_logical <- rep(FALSE, nparticles) # if there is stability issues
    accept_logical[finite_logical] <- (log(runif(nfinite)) < accept_ratio[finite_logical]) # if there is stability issues
    # accept_logical <- (log(runif(nparticles)) < accept_ratio)
    x[accept_logical, ] <- proposed_x[accept_logical, ]
    logtarget_x[accept_logical] <- logtarget_proposal[accept_logical]
    acceptprob <- sum(accept_logical) / nparticles
    
    return(list(x = x, logtarget_x = logtarget_x, acceptprob = acceptprob))
  }
  
  return(list(kernel = kernel))
}
