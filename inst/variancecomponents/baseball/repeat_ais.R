rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(baseball_artificial_logprior(x))
prior$gradlogdensity <- function(x) baseball_gradlogprior_artificial(x)
prior$rinit <- function(n) baseball_sample_artificial_prior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(baseball_logprior(x) + 
                                                  baseball_loglikelihood(x) -
                                                  baseball_artificial_logprior(x))
likelihood$gradlogdensity <- function(x) baseball_gradlogprior(x) + 
  baseball_gradloglikelihood(x) - 
  baseball_gradlogprior_artificial(x)

# SMC settings
nparticles <- 2^7
nsteps <- 50
lambda <- seq(0, 1, length.out = nsteps)^2
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.05
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

save(ess, log_normconst, file = "inst/variancecomponents/baseball/results/repeat_ais.RData", safe = F)

