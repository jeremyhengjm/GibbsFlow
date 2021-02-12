rm(list = ls())
library(GibbsFlow)
library(tictoc)
library(ggplot2)
library(ggthemes)
library(pracma)
library(gridExtra)
setmytheme()

# prior
prior <- list()
prior$logdensity <- function(x) as.numeric(gaussian_logprior(x))
prior$gradlogdensity <- function(x) gaussian_gradlogprior(x)
prior$rinit <- function(n) gaussian_sampleprior(n) 

# likelihood
likelihood <- list()
likelihood$logdensity <- function(x) as.numeric(gaussian_loglikelihood(x))
likelihood$gradlogdensity <- function(x) gaussian_gradloglikelihood(x)
compute_gibbsflow <- function(stepsize, lambda, lambda_next, derivative_lambda, x, logdensity) gaussian_gibbsflow(stepsize, lambda, derivative_lambda, x, logdensity)

# posterior
posterior <- list()
posterior$mean <- function(lambda) gaussian_posterior_mean(lambda)
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)
posterior$log_normconst <- function(lambda) gaussian_log_normconst(lambda)
posterior$logdensity <- function(x, lambda) gaussian_logposterior(x, lambda)
posterior$mean(1)
posterior$cov(1)

# numerical integration
dimension <- 4
exponent <- 2
bound <- 10
gibbs_velocity <- function(t, x){
  # cat(x, "\n")
  if (any(abs(x) > bound)){
    cat("Particle outside boundary")
    output <- matrix(0, nrow = 1, ncol = dimension)
  } else {
    output <- as.matrix(gaussian_gibbsvelocity(t, x, exponent))
  }
  return(output)
}

# plot stepsize taken by adaptive integrator for four trajectories
nrepeats <- 4
plot.df <- data.frame()
for (i in 1:nrepeats){
  # initialize
  initial_condition <- as.numeric(prior$rinit(1))
  
  # run numerical integrator ode45
  output <- ode45(f = gibbs_velocity, t0 = 0, tfinal = 1, y0 = initial_condition)
  # atol = 1e-8, hmax = 1e-3)
  
  # examine stepsize 
  stepsizes <- diff(output$t)
  min_step <- min(stepsizes)
  cat("Minimum stepsize taken by ode45:", min_step, "\n")
  
  # examine number of steps
  nsteps <- length(output$t)
  cat("Number of steps taken by ode45:", nsteps, "\n")
  
  # store stuff for plotting
  plot.df <- rbind(plot.df, data.frame(time = output$t, stepsize = c(0, stepsizes),
                                       trajectory = factor(rep(i, nsteps))))
}
gtrajectory <- ggplot(data = plot.df) + geom_line(aes(x = time, y = stepsize, colour = trajectory)) +
  scale_color_colorblind() 
gtrajectory
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_stepsizes_demo.eps", plot = gtrajectory,
       device = "eps", width = 6, height = 6)

# repeats to investigate number of time steps and minimum step size taken by adaptive integrator
nrepeats <- 2^10
min_stepsize <- rep(0, nrepeats)
total_nsteps <- rep(0, nrepeats)
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

}

# plotting 
boxplots <- list()

# total number of time steps plot
timesteps.df <- data.frame(timesteps = total_nsteps)
boxplots[[1]] <- ggplot(timesteps.df, aes(x = NULL, y = timesteps)) + geom_boxplot() + coord_flip() + scale_x_discrete(labels = c()) + 
  xlab("") + ylab("") + ggtitle("no. of time steps") + theme(plot.title = element_text(hjust = 0.5))

# minimum stepsize plot
default_stepsize <- 0.01
cat("Proportion of trajectories less than default stepsize:", mean(min_stepsize != default_stepsize), "\n")
stepsize.df <- data.frame(stepsize = min_stepsize[min_stepsize != default_stepsize])
boxplots[[2]] <- ggplot(stepsize.df, aes(x = NULL, y = stepsize)) + geom_boxplot() + coord_flip() + scale_x_discrete(labels = c()) + 
    xlab("") + ylab("") + ggtitle("minimum stepsize") + theme(plot.title = element_text(hjust = 0.5))

# combine boxplots
combine_boxplots <- arrangeGrob(boxplots[[1]], boxplots[[2]],
                                ncol = 1, nrow = 2)
gboxplots <- grid.arrange(combine_boxplots)
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_stiffness_demo.eps", plot = gboxplots,
       device = "eps", width = 6, height = 6)

