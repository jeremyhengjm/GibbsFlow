rm(list = ls())
library(GibbsFlow)
library(spatstat)
library(tictoc)
library(ggplot2)

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

# artificial prior
prior <- list()
prior$logdensity <- function(x) as.numeric(coxprocess_log_ep_proposal(x))
prior$gradlogdensity <- function(x) coxprocess_gradlog_ep_proposal(x)
prior$rinit <- function(n) coxprocess_sample_ep_proposal(n)

# artificial likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(coxprocess_logprior(x)) + 
  as.numeric(coxprocess_loglikelihood(x, data_counts)) - 
  as.numeric(coxprocess_log_ep_proposal(x))
likelihood$gradlogdensity <- function(x) coxprocess_gradlogprior(x) + 
  coxprocess_gradloglikelihood(x, data_counts) - 
  coxprocess_gradlog_ep_proposal(x)

# smc settings
nparticles <- 2^9
nsteps <- 80
lambda <- seq(0, 1, length.out = nsteps)^2
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.25
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 2

# run ais/smc
tic()
smc <- run_ais(prior, likelihood, nparticles, lambda, mcmc)
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

# acceptance probability
acceptprob_min.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[1, ])
acceptprob_max.df <- data.frame(time = 1:nsteps, acceptprob = smc$acceptprob[2, ])
ggplot() + geom_line(data = acceptprob_min.df, aes(x = time, y = acceptprob), colour = "blue") + 
  geom_line(data = acceptprob_max.df, aes(x = time, y = acceptprob), colour = "red") + ylim(c(0, 1))

