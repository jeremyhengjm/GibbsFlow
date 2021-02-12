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

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(coxprocess_logprior(x))
prior$gradlogdensity <- function(x) coxprocess_gradlogprior(x)
prior$rinit <- function(n) coxprocess_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(coxprocess_loglikelihood(x, data_counts))
likelihood$gradlogdensity <- function(x) coxprocess_gradloglikelihood(x, data_counts)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity, counts) coxprocess_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity, data_counts)

# smc settings
nparticles <- 2^9
nsteps <- 80
timegrid <- seq(0, 1, length.out = nsteps)
exponent <- 2
lambda <- timegrid^exponent
derivative_lambda <- exponent * timegrid^(exponent - 1)

# repeat gibbsflow sis/sisr
nrepeats <- 100
ess <- matrix(nrow = nrepeats, ncol = nsteps)
log_normconst <- rep(0, nrepeats)
for (irepeat in 1:nrepeats){
  cat("Repeat:", irepeat, "/", nrepeats, "\n")
  tic()
    smc <- run_gibbsflow_sis(prior, likelihood, nparticles, timegrid, lambda, derivative_lambda, compute_gibbsflow)
  toc()
  ess[irepeat, ] <- smc$ess * 100 / nparticles
  log_normconst[irepeat] <- smc$log_normconst[nsteps]
}

save(ess, log_normconst, file = "inst/coxprocess/results/repeat_gibbsflow_sis_400.RData", safe = F)
