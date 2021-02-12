rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
library(ggthemes)
library(deSolve)
setmytheme()

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(mixtureexample_logprior(x))
prior$gradlogdensity <- function(x) mixtureexample_gradlogprior(x) 
prior$rinit <- function(n) mixtureexample_sampleprior(n)

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(mixtureexample_loglikelihood(x))
likelihood$gradlogdensity <- function(x) mixtureexample_gradloglikelihood(x)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) mixtureexample_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity)

# posterior
posterior <- list()
posterior$mean <- function(component) mixtureexample_posterior_mean(component-1)
posterior$cov <- function(component) mixtureexample_posterior_cov(component-1)
posterior$log_normconst <- mixtureexample_log_normconst()
posterior$logdensity <- function(x) mixtureexample_logposterior(x)

# gibbs velocity
exponent <- 2
gibbs_velocity <- function(t, x, parms){
  output <- list(as.numeric(mixtureexample_gibbsvelocity(t, x, exponent)))
  return(output)
}

# demo compute gibbsflow
initial_condition <- as.numeric(prior$rinit(1))
times <- seq(0, 1, length.out = 11)
tic()
output_ode <- ode(y = initial_condition, times = times, func = gibbs_velocity)#, method = "radau")
toc()
output_ode

# repeat computation of gibbsflow
nparticles <- 2^14
dimension <- 2
smc <- list()
smc$xparticles <- matrix(nrow = nparticles, ncol = dimension)
times <- c(0, 1)
tic()
for (n in 1:nparticles){
  initial_condition <- as.numeric(prior$rinit(1))
  output_ode <- ode(y = initial_condition, times = times, func = gibbs_velocity)#, method = "radau")
  smc$xparticles[n, ] <- output_ode[2,2:(dimension+1)]
  cat("Particle: ", n, "/", nparticles, "\n")
  cat("Terminal position: ", smc$xparticles[n, ], "\n")
  if ( (n %% 1000) == 0 ){
    save(smc, file = "inst/mixtureexample/results/xi_6_rho_75.RData")
  }
}
toc()
like_xi <- 6
like_rho <- 0.75
save(like_xi, like_rho, smc, file = "inst/mixtureexample/results/xi_6_rho_75.RData")

# plot density function and trajectory
# load("inst/mixtureexample/results/xi_6_rho_0.RData")
# load("inst/mixtureexample/results/xi_6_rho_15.RData")
# load("inst/mixtureexample/results/xi_6_rho_30.RData")
# load("inst/mixtureexample/results/xi_6_rho_45.RData")
# load("inst/mixtureexample/results/xi_6_rho_60.RData")
# load("inst/mixtureexample/results/xi_6_rho_75.RData")
load("inst/mixtureexample/results/xi_6_rho_90.RData")

# plot proportion of particles in each mode
nparticles <- 2^14
cluster_particles <- kmeans(smc$xparticles, t(sapply(1:4, function(i) posterior$mean(i))))
cat("Empirical proportion:", table(cluster_particles$cluster) / nparticles)
cluster_particles$centers # check centers
t(sapply(1:4, function(i) posterior$mean(i))) # posterior means

# Pearson's Chi-squared test for uniformity 
chisq.test(x = cluster_particles$size, p = c(0.4, 0.1, 0.4, 0.1))

# plot terminal positions
xgrid <- seq(-10, 10, length.out = 200)
ygrid <- seq(-10, 10, length.out = 200)
grid <- expand.grid(xgrid = xgrid, ygrid = ygrid)
prior_values <- exp(prior$logdensity(as.matrix(grid)))
posterior_values <- exp(posterior$logdensity(as.matrix(grid)))
plot_prior <- cbind(grid, prior_values)
plot_posterior <- cbind(grid, posterior_values)
plot_particles <- data.frame(x = smc$xparticles[, 1], y = smc$xparticles[, 2], 
                             cluster = as.character(cluster_particles$cluster))
gcontour <- ggplot() + stat_contour(data = plot_prior, aes(x = xgrid, y = ygrid, z = prior_values), colour = "red") +
  stat_contour(data = plot_posterior, aes(x = xgrid, y = ygrid, z = posterior_values)) +
  # geom_point(data = plot_particles, aes(x = x, y = y, colour = cluster), size = 1, alpha = 0.4) + 
  geom_point(data = plot_particles, aes(x = x, y = y, colour = cluster), size = 1, alpha = 1) + 
  scale_colour_colorblind() + xlim(-6, 6) + ylim(-6, 6) + 
  coord_fixed(ratio = 1) + labs(x = expression(x[1]), y = expression(x[2]), title = bquote(rho == 0.90))
  
gcontour
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/arXiv_v2/gmixture_rho_90.eps", plot = gcontour,
       device = "eps", width = 6, height = 6)

