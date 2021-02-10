#' @rdname get_mala_kernel
#' @title Construct Metropolis adjusted Langevin algorithm kernel
#' @param logtarget evaluates log density of the target distribution  
#' @param gradlogtarget evaluates gradient of the log target density 
#' @param parameters list with key: 
#' \code{stepsize} specifies stepsize of discretization
#' @return list with key: \code{kernel}
#' @export
get_mala_kernel <- function(logtarget, gradlogtarget, parameters){
  stepsize <- parameters$stepsize
  stepsize_squared <- stepsize^2
  
  # One step of MALA
  kernel <- function(x, logtarget_x){
    # x: N x d matrix
    nparticles <- nrow(x)
    dimension <- ncol(x)
    
    # compute forward move
    forward_euler <- x + 0.5 * stepsize_squared * gradlogtarget(x)
    
    # sample proposal
    increment <- matrix(rnorm(nparticles * dimension), nrow = nparticles, ncol = dimension) # could be faster?
    proposal <- forward_euler + stepsize * increment
    
    # compute backward move
    backward_euler <- proposal + 0.5 * stepsize_squared * gradlogtarget(proposal)
    
    # compute acceptance probability
    logtarget_proposal <- logtarget(proposal) 
    accept_ratio <- logtarget_proposal - logtarget_x - 
      0.5 * rowSums((x - backward_euler)^2) / stepsize_squared + 
      0.5 * rowSums((proposal - forward_euler)^2) / stepsize_squared
    
    # accept or reject proposal
    accept_logical <- (log(runif(nparticles)) < accept_ratio)
    x[accept_logical, ] <- proposal[accept_logical, ]
    logtarget_x[accept_logical] <- logtarget_proposal[accept_logical]
    acceptprob <- sum(accept_logical) / nparticles
    
    return(list(x = x, logtarget_x = logtarget_x, acceptprob = acceptprob))
  }
  
  return(list(kernel = kernel))
}
