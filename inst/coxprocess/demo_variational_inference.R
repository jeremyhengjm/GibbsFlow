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

# prior distribution
parameter_sigmasq <- 1.91
parameter_invsigmasq <- 1 / parameter_sigmasq;
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
parameter_beta <- 1 / 33
parameter_area <- 1 / dimension

prior_mean <- rep(parameter_mu, dimension)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for (m in 1:dimension){
  for (n in 1:dimension){
    index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
    index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
    prior_cov[m,n] <- exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameter_beta) ) # excludes parameter_sigmasq factor
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))

# posterior distribution
posterior_logdensity <- function(x) as.numeric(coxprocess_logprior(x)) + 
  as.numeric(coxprocess_loglikelihood(x, data_counts))

# mean field variational bayes as describe in Teng, Nathoo & Johnson (2017)
uniroot_range <- c(parameter_mu - 3 * sqrt(parameter_sigmasq), parameter_mu + 3 * sqrt(parameter_sigmasq))
niterations <- 32
vb_mean <- prior_mean
vb_var <- rep(parameter_sigmasq, dimension)
vb_mean_old <- vb_mean
vb_var_old <- vb_var 

inv_hessians <- 1 / (parameter_area * exp(vb_mean) + parameter_invsigmasq * diag(prior_precision))
elbo_old <- sum(vb_mean * data_counts - parameter_area * exp(vb_mean + 0.5 * inv_hessians)) + 0.5 * sum(log(inv_hessians))

tic()
for (iteration in 1:niterations){
  cat("Iteration:", iteration, "\n")
  for (i in 1:dimension){
    # define f function in equation 8 
    f <- function(x){
      copy_vb_mean <- vb_mean
      copy_vb_mean[i] <- x
      return(data_counts[i] - parameter_area * exp(x) - parameter_invsigmasq * sum(prior_precision[i, ] * (copy_vb_mean - parameter_mu)))
    }
    
    # Laplace approximation 
    new_mean <- uniroot(f, uniroot_range, extendInt = "yes")$root
    vb_mean[i] <- new_mean
    
    new_var <- parameter_area * exp(new_mean) + parameter_invsigmasq * prior_precision[i, i] # negative Hessian 
    vb_var[i] <- 1 / new_var
 
  }
  
  cat("Change in means:", sum(abs(vb_mean - vb_mean_old)), "\n")
  cat("Change in variances", sum(abs(vb_var - vb_var_old)), "\n")
  vb_mean_old <- vb_mean
  vb_var_old <- vb_var
  
  # compute evidence lower bound
  inv_hessians <- 1 / (parameter_area * exp(vb_mean) + parameter_invsigmasq * diag(prior_precision))
  elbo <- sum(vb_mean * data_counts - parameter_area * exp(vb_mean + 0.5 * inv_hessians)) + 0.5 * sum(log(inv_hessians))
  cat("Change in ELBO:", abs(elbo - elbo_old), "\n")
  elbo_old <- elbo
  
}
toc()

save(vb_mean, vb_var, file = "inst/coxprocess/vi_proposal.RData")

# demo importance sampling from variational distribution
nparticles <- 32
samples <- mvnrnd(nparticles, vb_mean, diag(vb_var))
logproposal <- as.numeric(mvnpdf(samples, vb_mean, diag(vb_var)))
logweights <- posterior_logdensity(samples) - logproposal
maxlogweights <- max(logweights)
weights <- exp(logweights - maxlogweights)
normweights <- weights / sum(weights)
ess <- 1 / sum(normweights^2)
ess / nparticles * 100
log_normconst <- log(mean(weights)) + maxlogweights

# repeat importance sampling from variational distribution
nrepeats <- 100
nparticles <- 32
ess <- rep(0, nrepeats)
log_normconst <- rep(0, nrepeats)
tic()
for (irepeat in 1:nrepeats){
  samples <- mvnrnd(nparticles, vb_mean, diag(vb_var))
  logproposal <- as.numeric(mvnpdf(samples, vb_mean, diag(vb_var)))
  logweights <- posterior_logdensity(samples) - logproposal
  maxlogweights <- max(logweights)
  weights <- exp(logweights - maxlogweights)
  normweights <- weights / sum(weights)
  ess[irepeat] <- (1 / sum(normweights^2)) * 100 / nparticles
  log_normconst[irepeat] <- log(mean(weights)) + maxlogweights
}
toc()

cat("Average ESS%:", mean(ess), "\n")
cat("Variance of log normalizing constant estimator:", var(log_normconst), "\n")
