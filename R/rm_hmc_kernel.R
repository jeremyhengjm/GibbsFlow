#'@rdname get_rm_hmc_kernel
#'@title Construct Riemann manifold Hamiltonian Monte Carlo kernel for Cox process application
#'@param logtarget evaluates log density of the target distribution  
#'@param gradlogtarget evaluates gradient of the log target density 
#' @param parameters list with keys: 
#' \code{stepsize} specifies stepsize of leap-frog discretization, 
#' \code{nsteps} specifies number of leap-frog steps
#' \code{metric} specifies metric tensor
#' @return list with key: \code{kernel}
#'@export
get_rm_hmc_kernel <- function(logtarget, gradlogtarget, parameters){
  stepsize <- parameters$stepsize
  nsteps <- parameters$nsteps
  metric <- parameters$metric
  
  # leap frog integrator
  leapfrog <- function(x, v){
    # x: N x d matrix; v: N x d matrix
    v <- v + stepsize * gradlogtarget(x) / 2
    for (istep in 1:nsteps){
      x <- x + stepsize * (v %*% metric$inverse)
      if (istep != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    return(list(x = x, v = v))
  }
  
  # One step of RM-HMC
  kernel <- function(x, logtarget_x){
    # x: N x d matrix
    nparticles <- nrow(x)
    dimension <- ncol(x)
    
    # sample velocity
    v <- matrix(rnorm(nparticles * dimension), nrow = nparticles, ncol = dimension) %*% 
         t(metric$inverse_chol_inverse) 
    
    # run leap-frog integrator
    leapfrog_result <- leapfrog(x, v)
    proposed_x <- leapfrog_result$x
    proposed_v <- leapfrog_result$v
    
    # compute acceptance probability
    logtarget_proposal <- logtarget(proposed_x)
    accept_ratio <- logtarget_proposal - logtarget_x 
    accept_ratio <- accept_ratio + 0.5 * rowSums((v %*% metric$inverse) * v) -
                                   0.5 * rowSums((proposed_v %*% metric$inverse) * proposed_v) 
    accept_logical <- (log(runif(nparticles)) < accept_ratio)
    x[accept_logical, ] <- proposed_x[accept_logical, ]
    logtarget_x[accept_logical] <- logtarget_proposal[accept_logical]
    acceptprob <- sum(accept_logical) / nparticles
    
    return(list(x = x, logtarget_x = logtarget_x, acceptprob = acceptprob))
    
  }
  
  return(list(kernel = kernel))
}
