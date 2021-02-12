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

# define functions to compute gibbs flow (and optionally velocity)
exponent <- 2
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) baseball_gibbsflow(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity)
gibbsvelocity <- function(t, x) as.matrix(baseball_gibbsvelocity(t, x, exponent))
  
# smc settings
nparticles <- 2^7
nsteps <- 50
timegrid <- seq(0, 1, length.out = nsteps)
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# run sampler
tic()
smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, gibbsvelocity)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# normalizing constant plot
normconst.df <- data.frame(time = 1:nsteps, normconst = smc$log_normconst)
ggplot() + geom_line(data = normconst.df, aes(x = time, y = normconst), colour = "blue") +
  labs(x = "time", y = "log normalizing constant")

# norm of gibbs velocity
normvelocity.df <- data.frame(time = timegrid, 
                              lower = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.25)),
                              median = apply(smc$normvelocity, 2, median),
                              upper = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.75)))
gnormvelocity <- ggplot(normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 0.5) + 
  xlim(0, 1) + # scale_y_continuous(breaks = c(0, 40, 80, 120)) +
  xlab("time") + ylab("norm of Gibbs velocity") 
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_baseball_normvelocity_gfsis.pdf", plot = gnormvelocity,
       device = "pdf", width = 6, height = 6)

