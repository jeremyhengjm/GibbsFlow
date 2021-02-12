rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(pracma)
setmytheme()

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(gaussian_logprior(x))
prior$gradlogdensity <- function(x) gaussian_gradlogprior(x)
prior$rinit <- function(n) gaussian_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(gaussian_loglikelihood(x))
likelihood$gradlogdensity <- function(x) gaussian_gradloglikelihood(x)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) gaussian_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)
posterior$mean(1)
posterior$cov(1)

# numerical integration
dimension <- 2
exponent <- 2
bound <- 20
gibbs_velocity <- function(t, x){
  # cat(x, "\n")
  if (any(abs(x) > bound)){
    cat("Particle outside boundary")
    output <- matrix(0, nrow = 1, ncol = dimension)
  } else {
    output <- as.matrix(gaussian_gibbsvelocity(t, x, exponent))
  }
  return(output)
}

# repeats to investigate number of time steps and minimum step size taken by adaptive integrator
nrepeats <- 2^9
min_stepsize <- rep(0, nrepeats)
total_nsteps <- rep(0, nrepeats)
tic()
for (i in 1:nrepeats){
  cat("Repeat", i, "/", nrepeats, "\n")
  # initialize
  initial_condition <- as.numeric(prior$rinit(1))
  
  # run numerical integrator ode45
  output <- ode45(f = gibbs_velocity, t0 = 0, tfinal = 1, y0 = initial_condition)
  # atol = 1e-8, hmax = 1e-3)
  
  # examine stepsize 
  stepsizes <- diff(output$t)
  min_step <- min(stepsizes)
  min_stepsize[i] <- min_step
  cat("Minimum stepsize taken by ode45:", min_step, "\n")
  
  # examine number of steps
  nsteps <- length(output$t)
  total_nsteps[i] <- nsteps
  cat("Number of steps taken by ode45:", nsteps, "\n")
  
}
toc()
# save(min_stepsize, total_nsteps, file = "inst/gaussian/results/integrator_y0.RData", safe = F)
# save(min_stepsize, total_nsteps, file = "inst/gaussian/results/integrator_y10.RData", safe = F)
# save(min_stepsize, total_nsteps, file = "inst/gaussian/results/integrator_y20.RData", safe = F)
# save(min_stepsize, total_nsteps, file = "inst/gaussian/results/integrator_y30.RData", safe = F)
# save(min_stepsize, total_nsteps, file = "inst/gaussian/results/integrator_y40.RData", safe = F)
