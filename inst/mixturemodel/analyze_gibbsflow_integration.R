rm(list=ls())
library(ggplot2)
library(ggthemes)
library(gridExtra)
setmytheme()
load("inst/mixturemodel/gibbsflow_integration.RData")
nparticles <- length(total_nsteps)
# summarize total number of numerical integration steps
summary(total_nsteps)

# summarize minimum step size over unit time interval 
default_stepsize <- 0.01
non_default <- (min_stepsize != default_stepsize) # remove default step size
total_non_default <- sum(non_default)
cat("Proportion of non-default stepsizes:", total_non_default / nparticles, "\n")
summary(min_stepsize[non_default])

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
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_stiffness_demo.eps", plot = gboxplots,
       # device = "eps", width = 6, height = 6)
  
# plot norm of Gibbs velocity over time for GF-SIS particles 
rm(list=ls())
load("inst/mixturemodel/gibbsflow_trajectory_new.RData")
nsteps <- ncol(smc$normvelocity)
normvelocity.df <- data.frame(time = seq(0, 1, length.out = nsteps), 
                              lower = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.25)),
                              median = apply(smc$normvelocity, 2, median),
                              upper = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.75)))
gnormvelocity <- ggplot(normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 0.5) + 
  xlim(0, 1) + scale_y_continuous(breaks = c(0, 40, 80, 120)) +
  xlab("time") + ylab("norm of Gibbs velocity") 
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_normvelocity_gfsis.pdf", plot = gnormvelocity,
       device = "pdf", width = 6, height = 6)

# plot norm of Gibbs velocity over time for GF-AIS particles
rm(list=ls())
load("inst/mixturemodel/gibbsflow_with_hmc_new.RData")
nsteps <- ncol(smc$normvelocity)
normvelocity.df <- data.frame(time = seq(0, 1, length.out = nsteps), 
                              lower = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.25)),
                              median = apply(smc$normvelocity, 2, median),
                              upper = apply(smc$normvelocity, 2, function(x) quantile(x, probs = 0.75)))
gnormvelocity <- ggplot(normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 0.5) + 
  xlim(0, 1) + scale_y_continuous(breaks = c(0, 40, 80, 120)) +
  xlab("time") + ylab("norm of Gibbs velocity") 
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_normvelocity_gfais.pdf", plot = gnormvelocity,
       device = "pdf", width = 6, height = 6)
