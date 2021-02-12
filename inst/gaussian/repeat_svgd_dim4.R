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

# svgd settings (d = 1)
nfparticles <- 2^9
nlparticles <- 2^9
niterations <- 25
stepsize <- 0.1
bandwidth <- 25

# repeat SVGD
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = niterations)
log_normconst <- rep(0, nrepeats)
for (irepeat in 1:nrepeats){
  cat("SVGD repeat:", irepeat, "/", nrepeats, "\n")
  tic()
  svgd <- gaussian_stein_variational_importance_sampling(nfparticles, nlparticles, niterations, stepsize, bandwidth)
  toc()
  ess[irepeat, ] <- svgd$ess
  log_normconst[irepeat] <- svgd$log_normconst
}

save(ess, log_normconst, file = "inst/gaussian/results/repeat_svgd_dim4.RData", safe = F)
