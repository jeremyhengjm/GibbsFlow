rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(gaussian_logprior(x))
prior$gradlogdensity <- function(x) gaussian_gradlogprior(x)
prior$rinit <- function(n) gaussian_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(gaussian_loglikelihood(x))
likelihood$gradlogdensity <- function(x) gaussian_gradloglikelihood(x)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)

# SMC settings
nparticles <- 2^9
nsteps <- 100
lambda <- seq(0, 1, length.out = nsteps)^2
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.25
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 50

# repeat gibbsflow ais/smc
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
for (irepeat in 1:nrepeats){
  cat("AIS repeat:", irepeat, "/", nrepeats, "\n")
  tic()
  smc <- run_ais(prior, likelihood, nparticles, lambda, mcmc)
  toc()
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}
  
save(ess, log_normconst, file = "inst/gaussian/results/repeat_ais_dim4.RData", safe = F)