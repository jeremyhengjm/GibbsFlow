rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# problem specification
nobservations <- 100
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
load("inst/variancecomponents/simulated/dataset_102.RData")
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

# SMC settings
nparticles <- 2^7
nsteps <- 500
lambda <- seq(0, 1, length.out = nsteps)^2
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.0075
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 2

# repeat ais/smc
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
tic()
for (irepeat in 1:nrepeats){
  cat("Repeat:", irepeat, "/", nrepeats, "\n")
  smc <- run_ais(prior, likelihood, nparticles, lambda, mcmc)
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}
toc()

save(ess, log_normconst, file = "inst/variancecomponents/simulated/results/repeat_ais_102.RData", safe = F)

