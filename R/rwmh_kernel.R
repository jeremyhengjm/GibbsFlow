#' @rdname get_rwmh_kernel
#' @title Construct random walk Metropolis-Hastings kernel
#' @param logtarget evaluates log density of the target distribution  
#' @param parameters list with key: 
#' \code{stepsize} specifies stepsize of random walk
#' @return list with key: \code{kernel}
#' @export
get_rwmh_kernel <- function(logtarget, parameters){
  stepsize <- parameters$stepsize
  stepsize_squared <- stepsize^2
  
  # One step of RWMH
  kernel <- function(x, logtarget_x){
    # x: N x d matrix
    nparticles <- nrow(x)
    dimension <- ncol(x)
    
    # sample proposal
    increment <- matrix(rnorm(nparticles * dimension), nrow = nparticles, ncol = dimension) # could be faster?
    proposal <- x + stepsize * increment
    
    # compute acceptance probability
    logtarget_proposal <- logtarget(proposal) 
    accept_ratio <- logtarget_proposal - logtarget_x 
    
    # accept or reject proposal
    accept_logical <- (log(runif(nparticles)) < accept_ratio)
    x[accept_logical, ] <- proposal[accept_logical, ]
    logtarget_x[accept_logical] <- logtarget_proposal[accept_logical]
    acceptprob <- sum(accept_logical) / nparticles
    
    return(list(x = x, logtarget_x = logtarget_x, acceptprob = acceptprob))
  }
  
  return(list(kernel = kernel))
}
