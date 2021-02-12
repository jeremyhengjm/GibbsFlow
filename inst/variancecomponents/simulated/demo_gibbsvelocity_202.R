rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
library(pracma)

# problem specification
nobservations <- 200
nrepetitions <- 1
dimension <- nobservations + 2
sigma_e_sq <- 0.00434

# simulate dataset
# true_sigma_theta_sq <- 0.01
# true_mu <- 0
# true_theta <- true_mu + sqrt(true_sigma_theta_sq) * rnorm(nobservations)
# dataset <- matrix(nrow = nobservations, ncol = nrepetitions)
# for (i in 1:nobservations){
#   for (j in 1:nrepetitions){
#     dataset[i, j] <- true_theta[i] + sqrt(sigma_e_sq) * rnorm(1)
#   }
# }
# dataset_rowsums <- rowSums(dataset)

# load dataset 
load("inst/variancecomponents/simulated/dataset_202.RData")
dataset_rowsums <- rowSums(dataset)

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(varcomp_artificial_logprior(x, dimension))
prior$gradlogdensity <- function(x) varcomp_gradlogprior_artificial(x, dimension)
prior$rinit <- function(n) varcomp_sample_artificial_prior(n, dimension) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(varcomp_logprior(x, dimension) + 
                                                  varcomp_loglikelihood(x, dataset) -
                                                  varcomp_artificial_logprior(x, dimension))
likelihood$gradlogdensity <- function(x) varcomp_gradlogprior(x, dimension) + 
  varcomp_gradloglikelihood(x, dataset_rowsums, dimension, nrepetitions) - 
  varcomp_gradlogprior_artificial(x, dimension)

# gibbs velocity 
exponent <- 2
gibbs_velocity <- function(t, x){
  x <- matrix(x, nrow = 1)
  output <- as.matrix(varcomp_gibbsvelocity(t, x, exponent, dataset_rowsums, dimension))
  return(output)
}

# repeats
nrepeats <- 2^10
min_stepsize <- rep(0, nrepeats)
total_nsteps <- rep(0, nrepeats)
plot.df <- data.frame()
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
  
  # store stuff for plotting
  # plot.df <- rbind(plot.df, data.frame(time = output$t, stepsize = c(0, stepsizes),
  #                                      trajectory = factor(rep(i, nsteps))))
}
save(min_stepsize, total_nsteps, file = "inst/variancecomponents/simulated/results/gibbsflow_integration_202.RData", safe = F)

# g <- ggplot(data = plot.df) + geom_line(aes(x = time, y = stepsize, colour = trajectory)) +
#   scale_color_colorblind() 
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_stepsizes_y40.eps", plot = g,
#        device = "eps", width = 6, height = 6)

