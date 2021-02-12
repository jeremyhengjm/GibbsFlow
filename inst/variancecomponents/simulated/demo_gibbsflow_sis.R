rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# problem specification
nobservations <- 25
# nobservations <- 50
# nobservations <- 100
# nobservations <- 200
# nobservations <- 400
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
# load("inst/variancecomponents/simulated/dataset_52.RData")
# load("inst/variancecomponents/simulated/dataset_102.RData")
# load("inst/variancecomponents/simulated/dataset_202.RData")
# load("inst/variancecomponents/simulated/dataset_402.RData")
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

# define functions to compute gibbs flow (and optionally velocity)
exponent <- 2
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) varcomp_gibbsflow(stepsize, lambda, lambda_next, derivative_lambda, x, 
                                                                                                                 logdensity, dataset_rowsums, dimension)
gibbsvelocity <- function(t, x) as.matrix(varcomp_gibbsvelocity(t, x, exponent, dataset_rowsums, dimension))
  
# smc settings
nparticles <- 2^7
nsteps <- 125
# nsteps <- 250
# nsteps <- 500
# nsteps <- 750
# nsteps <- 1000
timegrid <- seq(0, 1, length.out = nsteps)
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# run sampler
tic()
# smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow)
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
save(normvelocity.df, file = "inst/variancecomponents/simulated/results/normvelocity_gfsis_27.RData")
# save(normvelocity.df, file = "inst/variancecomponents/simulated/results/normvelocity_gfsis_52.RData")
# save(normvelocity.df, file = "inst/variancecomponents/simulated/results/normvelocity_gfsis_102.RData")
# save(normvelocity.df, file = "inst/variancecomponents/simulated/results/normvelocity_gfsis_202.RData")
# save(normvelocity.df, file = "inst/variancecomponents/simulated/results/normvelocity_gfsis_402.RData")
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_baseball_normvelocity_gfais.pdf", plot = gnormvelocity,
# device = "pdf", width = 6, height = 6)


