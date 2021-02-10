#' @rdname construct_kernel
#' @title Constructs MCMC kernel
#' @param logtarget evaluates log density of the target distribution  
#' @param gradlogtarget evaluates gradient of the log target density 
#' @param mcmc list with keys: 
#' \code{choice} specifies type of MCMC method, 
#' \code{parameters} specifies algorithmic tuning parameters,
#' @return function which returns desired MCMC kernel
#' @export
construct_kernel <- function(logtarget, gradlogtarget, mcmc){
  if (mcmc$choice == "rwmh"){
    return(get_rwmh_kernel(logtarget, mcmc$parameters)$kernel)
  }
  
  if (mcmc$choice == "mala"){
    return(get_mala_kernel(logtarget, gradlogtarget, mcmc$parameters)$kernel)
  }
  
  if (mcmc$choice == "hmc"){
    return(get_hmc_kernel(logtarget, gradlogtarget, mcmc$parameters)$kernel)
  }
  
  if (mcmc$choice == "rmhmc"){
    return(get_rm_hmc_kernel(logtarget, gradlogtarget, mcmc$parameters)$kernel) # only for Cox process model
  }
  
}