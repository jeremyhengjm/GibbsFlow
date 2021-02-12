rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
setmytheme()

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(banana_logprior(x))
prior$gradlogdensity <- function(x) banana_gradlogprior(x) 
prior$rinit <- function(n) banana_sampleprior(n)

# likelihood
parameter_alpha <- 5
parameter_beta <- 10
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(banana_loglikelihood(x, parameter_alpha, parameter_beta))
likelihood$gradlogdensity <- function(x) banana_gradloglikelihood(x, parameter_alpha, parameter_beta)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) banana_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity, parameter_alpha, parameter_beta)

# smc settings
nparticles <- 2^8
nsteps <- 200
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# mcmc settings
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.05
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 1

# run sampler
set.seed(7)
tic()
smc <- run_gibbsflow_ais(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, mcmc)
# smc <- run_gibbsflow_smc(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, mcmc)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# normalizing constant plot
normconst.df <- data.frame(time = 1:nsteps, normconst = smc$log_normconst)
ggplot() + geom_line(data = normconst.df, aes(x = time, y = normconst), colour = "blue") + 
  labs(x = "time", y = "log normalizing constant") 

# acceptance probability
acceptprob_min.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[1, ])
acceptprob_max.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[2, ])
ggplot() + geom_line(data = acceptprob_min.df, aes(x = time, y = acceptprob), colour = "blue") + 
  geom_line(data = acceptprob_max.df, aes(x = time, y = acceptprob), colour = "red") + ylim(c(0, 1))

# plot density function and trajectory
posterior <- list()
posterior$logdensity <- function(x) banana_logprior(x) + banana_loglikelihood(x, parameter_alpha, parameter_beta) - smc$log_normconst[nsteps]
xgrid <- seq(-10, 10, length.out = 200)
ygrid <- seq(-10, 10, length.out = 200)
grid <- expand.grid(xgrid = xgrid, ygrid = ygrid)
prior_values <- exp(prior$logdensity(as.matrix(grid)))
posterior_values <- exp(posterior$logdensity(as.matrix(grid)))
plot_prior <- cbind(grid, prior_values)
plot_posterior <- cbind(grid, posterior_values)
plot_particles <- data.frame(x = smc$xparticles[, 1], y = smc$xparticles[, 2])
gcontour <- ggplot() + stat_contour(data = plot_prior, aes(x = xgrid, y = ygrid, z = prior_values), colour = "red") +
  stat_contour(data = plot_posterior, aes(x = xgrid, y = ygrid, z = posterior_values)) +
  geom_point(data = plot_particles, aes(x = x, y = y), colour = "black", size = 1) + 
  coord_fixed(ratio = 1) + labs(x = expression(x[1]), y = expression(x[2])) + 
  xlim(-3, 3) + ylim(-3, 6)
gcontour
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/banana_gibbsflow_ais.eps", plot = gcontour,
       device = "eps", width = 6, height = 6)
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/banana_gibbsflow_smc.eps", plot = gcontour,
       # device = "eps", width = 6, height = 6)

