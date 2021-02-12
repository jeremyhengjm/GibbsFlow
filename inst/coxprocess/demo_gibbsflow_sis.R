rm(list = ls())
library(GibbsFlow)
library(spatstat)
library(tictoc)
library(ggplot2)
setmytheme()

# load pine saplings dataset
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
plot(x = data_x, y = data_y, type = "p")

ngrid <- 20
grid <- seq(from = 0, to = 1, length.out = ngrid+1)
dimension <- ngrid^2
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(coxprocess_logprior(x))
prior$gradlogdensity <- function(x) coxprocess_gradlogprior(x)
prior$rinit <- function(n) coxprocess_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(coxprocess_loglikelihood(x, data_counts))
likelihood$gradlogdensity <- function(x) coxprocess_gradloglikelihood(x, data_counts)

# define functions to compute gibbs flow (and optionally velocity)
exponent <- 2
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) coxprocess_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity, data_counts)
gibbsvelocity <- function(t, x) as.matrix(coxprocess_gibbsvelocity(t, x, exponent, data_counts))

# smc settings
nparticles <- 2^9
# nsteps <- 40 # d = 10 x 10
# nsteps <- 60 # d = 15 x 15
nsteps <- 80 # d = 20 x 20
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# run sis/sisr
tic()
# smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow)
smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, gibbsvelocity)
toc()

# ess plot
ess.df <- data.frame(time = 1:nsteps, ess = smc$ess * 100 / nparticles)
ggplot(ess.df, aes(x = time, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))
smc$log_normconst[nsteps]

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
# save(normvelocity.df, file = "inst/coxprocess/results/normvelocity_gfsis_100.RData")
# save(normvelocity.df, file = "inst/coxprocess/results/normvelocity_gfsis_225.RData")
save(normvelocity.df, file = "inst/coxprocess/results/normvelocity_gfsis_400.RData")
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_normvelocity_gfsis_400.pdf", plot = gnormvelocity,
       # device = "pdf", width = 6, height = 6)

