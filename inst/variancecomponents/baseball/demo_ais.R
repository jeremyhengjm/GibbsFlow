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

# run SMC
tic("SMC runtime")
smc <- run_ais(prior, likelihood, nparticles, lambda, mcmc)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# normalizing constant plot
normconst.df <- data.frame(time = 1:nsteps, normconst = smc$log_normconst)
ggplot() + geom_line(data = normconst.df, aes(x = time, y = normconst), colour = "black") + 
  labs(x = "time", y = "log normalizing constant")

# acceptance probability
acceptprob_min.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[1, ])
acceptprob_max.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[2, ])
  ggplot() + geom_line(data = acceptprob_min.df, aes(x = time, y = acceptprob), colour = "blue") + 
  geom_line(data = acceptprob_max.df, aes(x = time, y = acceptprob), colour = "red") + ylim(c(0, 1))

  