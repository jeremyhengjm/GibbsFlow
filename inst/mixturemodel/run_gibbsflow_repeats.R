rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
library(deSolve)
library(gridExtra)
setmytheme()

# problem specification
dimension <- 4
parameter_sigma <- 0.55

# generate observations
set.seed(17)
true_means <- c(-3, 0, 3, 6)
mixtureweights <- rep(0.25, 4)
mixturelogweights <- log(mixtureweights)
nobservations <- 100
observations <- rep(0, nobservations)
for (i in 1:nobservations){
  observations[i] <- true_means[floor((i-1)/25) + 1] + parameter_sigma * rnorm(1)
}

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(mixturemodel_logprior(x))
prior$gradlogdensity <- function(x) mixturemodel_gradlogprior(x)
prior$rinit <- function(n) mixturemodel_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(mixturemodel_loglikelihood(x, observations, mixturelogweights))
likelihood$gradlogdensity <- function(x) mixturemodel_gradloglikelihood(x, observations, mixturelogweights)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity, obs) 
  mixturemodel_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity, observations = observations, mixturelogweights)

# gibbs velocity 
exponent <- 2
gibbs_velocity <- function(t, x, parms){
  output <- list(as.numeric(mixturemodel_gibbsvelocity(t, x, exponent, observations, mixturelogweights)))
  return(output)
}

# demo compute gibbsflow
initial_condition <- as.numeric(prior$rinit(1))
times <- seq(0, 1, length.out = 11)
tic()
output_ode <- ode(y = initial_condition, times = times, func = gibbs_velocity)#, method = "radau")
toc()
output_ode

# repeat computation of gibbsflow
nparticles <- 2^10
nsteps <- 400
times <- seq(0, 1, length.out = nsteps)
# times <- c(0, 1)
smc <- list()
smc$xtrajectory <- array(dim = c(nparticles, dimension, nsteps))
smc$normvelocity <- matrix(0, nparticles, nsteps)
tic()
for (n in 1:nparticles){
  # initialize
  initial_condition <- as.numeric(prior$rinit(1))
  
  # run numerical integrator
  output_ode <- ode(y = initial_condition, times = times, func = gibbs_velocity)#, method = "radau")
  smc$xtrajectory[n, , ] <- t(output_ode[,2:(dimension+1)])
  
  # store norm of velocity 
  for (istep in 1:nsteps){
    current_time <- times[istep]
    smc$normvelocity[n, istep] <- sqrt(sum(as.numeric(gibbs_velocity(current_time, output_ode[istep, 2:(dimension+1)], parms = NULL)[[1]])^2))
  }
  
  # print and save
  cat("Particle: ", n, "/", nparticles, "\n")
  cat("Terminal position: ", smc$xtrajectory[n, , nsteps], "\n")
  if ( (n %% 100) == 0 ){
    save(smc, file = "inst/mixturemodel/gibbsflow_trajectory_new.RData")
  }
  
}
toc()
save(smc, file = "inst/mixturemodel/gibbsflow_trajectory_new.RData")

# repeat computation of gibbsflow with HMC
nparticles <- 2^10
nsteps <- 400
times <- seq(0, 1, length.out = nsteps)
lambda <- times^exponent
smc <- list()
smc$xtrajectory <- array(dim = c(nparticles, dimension, nsteps))
smc$normvelocity <- matrix(0, nparticles, nsteps)
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.1
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 1

tic()
for (n in 1:nparticles){
  # initialize 
  initial_position <- as.numeric(prior$rinit(1))
  smc$xtrajectory[n, , 1] <- initial_position
  
  # store norm of velocity 
  smc$normvelocity[n, 1] <- sqrt(sum(as.numeric(gibbs_velocity(0, initial_position, parms = NULL)[[1]])^2))
  
  for (istep in 2:nsteps){
    # Gibbs flow
    time_interval <- times[(istep-1):istep]
    output_ode <- ode(y = initial_position, times = time_interval, func = gibbs_velocity)#, method = "radau")
    xparticle <- matrix(output_ode[2, 2:(dimension+1)], nrow = 1)
    
    # store norm of velocity 
    current_time <- times[istep]
    smc$normvelocity[n, istep] <- sqrt(sum(as.numeric(gibbs_velocity(current_time, xparticle, parms = NULL)[[1]])^2))
    
    # MCMC
    current_logtarget <- function(x) prior$logdensity(x) + lambda[istep] * likelihood$logdensity(x)
    current_gradlogtarget <- function(x) prior$gradlogdensity(x) + lambda[istep] * likelihood$gradlogdensity(x)
    transition_kernel <- construct_kernel(current_logtarget, current_gradlogtarget, mcmc)
    
    current_logdensity <- current_logtarget(xparticle)
    for (imove in 1:mcmc$nmoves){
      mcmc_output <- transition_kernel(xparticle, current_logdensity)
      xparticle <- mcmc_output$x
      current_logdensity <- mcmc_output$logtarget_x
    }
    
    initial_position <- xparticle 
    smc$xtrajectory[n, , istep] <- xparticle
    
  }
  
  cat("Particle: ", n, "/", nparticles, "\n")
  cat("Terminal position: ", xparticle, "\n")
  if ( (n %% 100) == 0 ){
    save(smc, file = "inst/mixturemodel/gibbsflow_with_hmc_new.RData")
  }
  
}
toc()
save(smc, file = "inst/mixturemodel/gibbsflow_with_hmc_new.RData")
