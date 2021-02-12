rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
library(ggthemes)
library(pracma)
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
gibbs_velocity <- function(t, x){
  output <- as.matrix(mixturemodel_gibbsvelocity(t, x, exponent, observations, mixturelogweights))
  return(output)
}

# plot stepsize taken by adaptive integrator for four trajectories
nrepeats <- 4
plot.df <- data.frame()
for (i in 1:nrepeats){
  # initialize
  initial_condition <- as.numeric(prior$rinit(1))
  
  # run numerical integrator ode45
  output <- ode45(f = gibbs_velocity, t0 = 0, tfinal = 1, y0 = initial_condition)
  # atol = 1e-8, hmax = 1e-3)
  
  # examine stepsize 
  stepsizes <- diff(output$t)
  min_step <- min(stepsizes)
  cat("Minimum stepsize taken by ode45:", min_step, "\n")
  
  # examine number of steps
  nsteps <- length(output$t)
  cat("Number of steps taken by ode45:", nsteps, "\n")
  
  # print terminal position
  cat(tail(output$y,1), "\n")
  
  # store stuff for plotting
  plot.df <- rbind(plot.df, data.frame(time = output$t, stepsize = c(0, stepsizes),
                                       trajectory = factor(rep(i, nsteps))))
}

gtrajectory <- ggplot(data = plot.df) + geom_line(aes(x = time, y = stepsize, colour = trajectory)) +
  scale_color_colorblind() 
gtrajectory
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_stepsizes.eps", plot = gtrajectory,
       device = "eps", width = 6, height = 6)

# repeats to investigate number of time steps and minimum step size taken by adaptive integrator
nrepeats <- 2^10
min_stepsize <- rep(0, nrepeats)
total_nsteps <- rep(0, nrepeats)
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
  
  # print terminal position
  cat(tail(output$y,1), "\n")
  
}

save(min_stepsize, total_nsteps, file = "inst/mixturemodel/gibbsflow_integration.RData", safe = F)
