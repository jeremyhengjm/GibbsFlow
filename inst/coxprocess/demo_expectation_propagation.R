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
    prior_cov[m,n] <- parameter_sigmasq * exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameter_beta) )
    
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))
prior_meanprecision <- prior_mean %*% prior_precision

# posterior distribution
posterior_logdensity <- function(x) as.numeric(coxprocess_logprior(x)) + 
  as.numeric(coxprocess_loglikelihood(x, data_counts))

# expectation propagation to compute Gaussian approximation of posterior
niterations <- 20
ep_mean <- prior_mean
ep_invsigmasq <- rep(parameter_invsigmasq, dimension)
ep_mean_old <- ep_mean
ep_invsigmasq_old <- ep_invsigmasq

totaltime <- 0
for (iteration in 1:niterations){
  tic()
  # update site i 
  for (i in 1:dimension){
    ep_invsigmasq_copy <- ep_invsigmasq
    ep_invsigmasq_copy[i] <- 0
    Q_lessi <- prior_precision + diag(ep_invsigmasq_copy)
    r_lessi <- prior_meanprecision + ep_mean * ep_invsigmasq_copy
    
    # define hybrid density and its gradient
    negative_logh <- function(x){
      output <- -0.5 * sum((x %*% Q_lessi) * x) + sum(x * r_lessi) + x[i] * data_counts[i] - parameter_area * exp(x[i])
      return(-output)
    }
    negative_gradlogh <- function(x){
      output <- - x %*% Q_lessi + r_lessi 
      output[i] <- output[i] + data_counts[i] - parameter_area * exp(x[i])
      return(-output)
    }
    
    # use Laplace approximation to compute mean and variance of hybrid distribution
    optim_output <- optim(par = r_lessi %*% solve(Q_lessi),
    # optim_output <- optim(par = r_lessi,
          fn = negative_logh,
          gr = negative_gradlogh, 
          method = "BFGS")
    mu_h <- as.numeric(optim_output$par)
    
    cat("Iteration:", iteration, "Dimension:", i,
    "Optim convergence", optim_output$convergence, "\n")
    invsigma_h <- Q_lessi 
    invsigma_h[i, i] <- invsigma_h[i, i] + parameter_area * exp(mu_h[i])
    
    # obtain natural parameters of site i
    Q_i <- invsigma_h - Q_lessi
    r_i <- mu_h %*% invsigma_h - r_lessi
    
    # reparameterize into mean and precision
    ep_invsigmasq[i] <- Q_i[i, i]
    ep_mean[i] <- r_i[i] / ep_invsigmasq[i]
    
  }
  
  cat("Iteration:", iteration, "Change in means:", sum(abs(ep_mean - ep_mean_old)), "\n")
  cat("Iteration:", iteration, "Change in variances", sum(abs(ep_invsigmasq - ep_invsigmasq_old)), "\n")
  ep_mean_old <- ep_mean
  ep_invsigmasq_old <- ep_invsigmasq
  
  timer <- toc()
  totaltime <- totaltime + (timer$toc - timer$tic)
  cat("Iteration", iteration, "Time elapsed", totaltime, "\n")
    
}

# compute Gaussian approximation 
ep_precision <- prior_precision + diag(ep_invsigmasq) # Q matrix
ep_r <- prior_meanprecision + ep_mean * ep_invsigmasq
ep_cov <- solve(ep_precision)
ep_meanvec <- ep_r %*% ep_cov

save(ep_meanvec, ep_cov, file = "inst/coxprocess/ep_proposal.RData")

# demo importance sampling from variational distribution
tic()
nparticles <- 100
samples <- mvnrnd(nparticles, ep_meanvec, ep_cov)
logproposal <- as.numeric(mvnpdf(samples, ep_meanvec, ep_cov))
logweights <- posterior_logdensity(samples) - logproposal
maxlogweights <- max(logweights)
weights <- exp(logweights - maxlogweights)
normweights <- weights / sum(weights)
ess <- 1 / sum(normweights^2)
ess / nparticles * 100
log_normconst <- log(mean(weights)) + maxlogweights
toc()

# repeat importance sampling from variational distribution
nrepeats <- 100
nparticles <- 32
ess <- rep(0, nrepeats)
log_normconst <- rep(0, nrepeats)
tic()
for (irepeat in 1:nrepeats){
  samples <- mvnrnd(nparticles, ep_meanvec, ep_cov)
  logproposal <- as.numeric(mvnpdf(samples, ep_meanvec, ep_cov))
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
