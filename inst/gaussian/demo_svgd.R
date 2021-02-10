rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(gaussian_logprior(x))
prior$gradlogdensity <- function(x) gaussian_gradlogprior(x)
prior$rinit <- function(n) gaussian_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(gaussian_loglikelihood(x))
likelihood$gradlogdensity <- function(x) gaussian_gradloglikelihood(x)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)

# svgd settings (d = 1)
nfparticles <- 2^9
nlparticles <- 2^9
niterations <- 45
stepsize <- 0.1
bandwidth <- 25

# run SVGD
tic()
svgd <- gaussian_stein_variational_importance_sampling(nfparticles, nlparticles, niterations, stepsize, bandwidth)
toc()

# ess plot
ess.df <- data.frame(iteration = 1:niterations, ess = svgd$ess)
ggplot(ess.df, aes(x = iteration, y = ess)) + geom_line() + 
  labs(x = "time", y = "ESS%") + ylim(c(0, 100))

# compare normalizing constant 
cat("True log-normalizing constant:", posterior$log_normconst(1), "\n")
cat("SVGD estimate:", svgd$log_normconst, "\n")

# histogram
hist(svgd$xparticles, probability = TRUE)
range <- seq(min(svgd$xparticles), max(svgd$xparticles), length.out = 100)
lines(x = range, y = dnorm(range, posterior$mean(1), sqrt(posterior$cov(1)[1, 1])), col = "red")

# # plot density function and trajectory
# xgrid <- seq(-3, 7, length.out = 100)
# ygrid <- seq(-3, 7, length.out = 100)
# grid <- expand.grid(xgrid = xgrid, ygrid = ygrid)
# prior_values <- rep(0, nrow(grid))
# posterior_values <- rep(0, nrow(grid))
# for (igrid in 1:nrow(grid)){
#   prior_values[igrid] <- exp(prior$logdensity(as.matrix(grid[igrid, ], nrow = 1)))
#   posterior_values[igrid] <- exp(posterior$logdensity(as.matrix(grid[igrid, ], nrow = 1), 1))
# }
# plot_prior <- cbind(grid, prior_values)
# plot_posterior <- cbind(grid, posterior_values)
# plot_particles <- data.frame(x = svgd$xparticles[, 1], y = svgd$xparticles[, 2])
# ggplot() + stat_contour(data = plot_prior, aes(x = xgrid, y = ygrid, z = prior_values), colour = "red") +
#   stat_contour(data = plot_posterior, aes(x = xgrid, y = ygrid, z = posterior_values)) +
#   geom_point(data = plot_particles, aes(x = x, y = y), colour = "black")
