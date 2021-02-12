rm(list=ls())
library(ggplot2)
library(ggthemes)
setmytheme()

nrepeats <- 2^10
default_stepsize <- 0.01
exclude_default_stepsize <- TRUE

# dimension 10 x 10
load("inst/coxprocess/results/gibbsflow_integration_100.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- data.frame(minstepsize = min_stepsize[non_default], 
                                dimension = factor(rep(100, total_non_default)))
  cat("Dimension:", 100, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- data.frame(minstepsize = min_stepsize, 
                                dimension = factor(rep(100, nrepeats)))
}
total.nsteps.df <- data.frame(timesteps = total_nsteps, 
                              dimension = factor(rep(100, nrepeats)))

# dimension 15 x 15
load("inst/coxprocess/results/gibbsflow_integration_225.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(225, total_non_default))))
  cat("Dimension:", 225, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(225, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(225, nrepeats))))

# dimension 20 x 20
load("inst/coxprocess/results/gibbsflow_integration_400.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(400, total_non_default))))
  cat("Dimension:", 400, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(400, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(400, nrepeats))))

# plot minimum step size against observation
gminstepsize <- ggplot(data = min.stepsize.df, aes(x = dimension, y = minstepsize)) + geom_boxplot() + scale_y_log10() + 
  xlab("dimension") + ylab("minimum stepsize")
gminstepsize
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_minstepsize_dim.eps", plot = gminstepsize,
       device = "eps", width = 6, height = 6)

# plot no. of time steps against observation
gtimesteps <- ggplot(data = total.nsteps.df, aes(x = dimension, y = timesteps)) + geom_boxplot() + 
  xlab("dimension") + ylab("time steps")
gtimesteps
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_timesteps_dim.eps", plot = gtimesteps,
       device = "eps", width = 6, height = 6)

