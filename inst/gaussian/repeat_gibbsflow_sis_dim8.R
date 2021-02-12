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
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) gaussian_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)

# smc settings
nparticles <- 2^9
nsteps <- 100
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# repeat gibbsflow sis
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
for (irepeat in 1:nrepeats){
  cat("GF-SIS repeat:", irepeat, "/", nrepeats, "\n")
  tic()
  smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow)
  toc()
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}

save(ess, log_normconst, file = "inst/gaussian/results/repeat_gibbsflow_sis_dim8.RData", safe = F)