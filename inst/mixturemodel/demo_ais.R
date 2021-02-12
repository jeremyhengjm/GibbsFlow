rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# problem specification
dimension <- 4
parameter_sigma <- 0.55

# generate observations
set.seed(17)
true_means <- c(-3, 0, 3, 6)
mixtureweights <- rep(0.25, 4)
mixturelogweights <- log(mixtureweights)

nobservations <- 100
observations <- rep(0, nobservations)
for (i in 1:nobservations){
  observations[i] <- true_means[floor((i-1)/25) + 1] + parameter_sigma * rnorm(1)
}

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(mixturemodel_logprior(x))
prior$gradlogdensity <- function(x) mixturemodel_gradlogprior(x)
prior$rinit <- function(n) mixturemodel_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(mixturemodel_loglikelihood(x, observations, mixturelogweights))
likelihood$gradlogdensity <- function(x) mixturemodel_gradloglikelihood(x, observations, mixturelogweights)

# smc settings
nparticles <- 2^14
nsteps <- 400
lambda <- seq(0, 1, length.out = nsteps)^2
mcmc <- list()
mcmc$choice <- "hmc"
mcmc$parameters$stepsize <- 0.1
mcmc$parameters$nsteps <- 10
mcmc$nmoves <- 10

# run ais/smc
tic()
smc_reference <- run_smc(prior, likelihood, nparticles, lambda, mcmc)
toc()

save(smc_reference, file = "inst/mixturemodel/smc_reference.RData", safe = F)
smc <- smc_reference

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

# plot particles
kde.df <- data.frame(x = smc$xparticles[, 1], y = smc$xparticles[, 2])
g <- ggplot(data = kde.df, aes(x = x, y = y)) + 
     geom_point(alpha = 0.5) + 
     geom_density_2d() + 
     xlim(-5, 8) + 
     ylim(-5, 8)
g

# plot proportion of particles in each mode
permutations <- function( x, prefix = c() )
{
  if(length(x) == 0 ) return(prefix)
  do.call(rbind, sapply(1:length(x), FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
}
cluster_particles <- kmeans(smc$xparticles, permutations(true_means))
# cat("Empirical proportion:", table(cluster_particles$cluster) / nparticles)
# cluster_particles$centers # check centers
barplot.df <- data.frame(index = cluster_particles$cluster)
g <- ggplot(barplot.df, aes(x = as.factor(index))) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  xlab("Mode") + 
  ylab("Proportion of particles")
g
