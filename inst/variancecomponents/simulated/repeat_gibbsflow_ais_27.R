rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# problem specification
nobservations <- 25
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
load("inst/variancecomponents/simulated/dataset_27.RData")
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
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) varcomp_gibbsflow(stepsize, lambda, lambda_next, derivative_lambda, x, 
                                                                                                                 logdensity, dataset_rowsums, dimension)

# smc settings
nparticles <- 2^7
nsteps <- 125
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# mcmc settings
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.0125
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 1

# repeat gibbsflow ais/smc
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
tic()
for (irepeat in 1:nrepeats){
  cat("Repeat:", irepeat, "/", nrepeats, "\n")
  smc <- run_gibbsflow_ais(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, mcmc)
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}
toc()

save(ess, log_normconst, file = "inst/variancecomponents/simulated/results/repeat_gibbsflow_ais_27.RData", safe = F)
