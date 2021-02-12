rm(list = ls())
library(GibbsFlow)
library(spatstat)
library(tictoc)

# load pine saplings dataset
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
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
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) coxprocess_gibbsflow_ep_proposal(stepsize, lambda, derivative_lambda, x, logdensity, data_counts)

# smc settings
nparticles <- 2^9
nsteps <- 80
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# mcmc settings
mcmc <- list()
mcmc$choice <- "rmhmc"
mcmc$parameters$stepsize <- 0.25
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 1

# compute metric tensor as in Girolami and Calderhead 2011 (Section 9)
parameter_sigmasq <- 1.91
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
parameter_beta <- 1 / 33
parameter_area <- 1 / dimension
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for (m in 1:dimension){
  for (n in 1:dimension){
    index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
    index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
    prior_cov[m,n] <- parameter_sigmasq * exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameter_beta) )
  }
}
prior_precision <- solve(prior_cov)
metric_tensor <- prior_precision
diag(metric_tensor) <- parameter_area * exp(parameter_mu + 0.5 * diag(prior_cov)) + diag(prior_precision)
mcmc$parameters$metric$inverse  <- solve(metric_tensor)
metric_chol_inverse <- t(chol(mcmc$parameters$metric$inverse))
mcmc$parameters$metric$inverse_chol_inverse <- solve(t(metric_chol_inverse))

# free up memory
prior_cov <- NULL
prior_precision <- NULL
metric_tensor <- NULL
metric_chol_inverse <- NULL

# repeat gibbsflow ais/smc
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
for (irepeat in 1:nrepeats){
  cat("Repeat:", irepeat, "/", nrepeats, "\n")
  tic()
    smc <- run_gibbsflow_ais(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow, mcmc)
  toc()
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}

save(ess, log_normconst, file = "inst/coxprocess/results/repeat_gibbsflow_ais_ep_400.RData", safe = F)
