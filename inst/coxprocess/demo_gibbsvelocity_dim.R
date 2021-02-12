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

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(coxprocess_logprior(x))
prior$gradlogdensity <- function(x) coxprocess_gradlogprior(x)
prior$rinit <- function(n) coxprocess_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(coxprocess_loglikelihood(x, data_counts))
likelihood$gradlogdensity <- function(x) coxprocess_gradloglikelihood(x, data_counts)

# gibbs velocity 
exponent <- 2
gibbs_velocity <- function(t, x){
  x <- matrix(x, nrow = 1)
  output <- as.matrix(coxprocess_gibbsvelocity(t, x, exponent, data_counts))
  return(output)
}

# repeats
nrepeats <- 2^10
min_stepsize <- rep(0, nrepeats)
total_nsteps <- rep(0, nrepeats)
plot.df <- data.frame()
for (i in 1:nrepeats){
  cat("Repeat", i, "/", nrepeats, "\n")
  # initialize
  initial_condition <- as.numeric(prior$rinit(1))
  
  # run numerical integrator ode45
  output <- ode45(f = gibbs_velocity, t0 = 0, tfinal = 1, y0 = initial_condition)
  # atol = 1e-8, hmax = 1e-3)
  
  # examine stepsize 
  stepsizes <- diff(output$t)
  min_step <- min(stepsizes)
  min_stepsize[i] <- min_step
  cat("Minimum stepsize taken by ode45:", min_step, "\n")
  
  # examine number of steps
  nsteps <- length(output$t)
  total_nsteps[i] <- nsteps
  cat("Number of steps taken by ode45:", nsteps, "\n")
  
  # store stuff for plotting
  # plot.df <- rbind(plot.df, data.frame(time = output$t, stepsize = c(0, stepsizes),
  #                                      trajectory = factor(rep(i, nsteps))))
}
# save(min_stepsize, total_nsteps, file = "inst/coxprocess/results/gibbsflow_integration_100.RData", safe = F)
# save(min_stepsize, total_nsteps, file = "inst/coxprocess/results/gibbsflow_integration_225.RData", safe = F)
save(min_stepsize, total_nsteps, file = "inst/coxprocess/results/gibbsflow_integration__400.RData", safe = F)

# g <- ggplot(data = plot.df) + geom_line(aes(x = time, y = stepsize, colour = trajectory)) +
#   scale_color_colorblind() # + ggtitle("quadratic schedule")
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_stepsizes_y40.eps", plot = g,
#        device = "eps", width = 6, height = 6)





