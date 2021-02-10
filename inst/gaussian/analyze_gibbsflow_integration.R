rm(list=ls())
library(ggplot2)
library(ggthemes)
setmytheme()

nrepeats <- 2^9
default_stepsize <- 0.01
exclude_default_stepsize <- TRUE

load("inst/gaussian/results/integrator_y0.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- data.frame(minstepsize = min_stepsize[non_default], 
                                observation = factor(rep(0, total_non_default)))
  cat("Proportion of non-default stepsize for xi =", 0, "is", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- data.frame(minstepsize = min_stepsize, 
                                observation = factor(rep(0, nrepeats)))
}
total.nsteps.df <- data.frame(timesteps = total_nsteps, 
                              observation = factor(rep(0, nrepeats)))

load("inst/gaussian/results/integrator_y10.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       observation = factor(rep(10, total_non_default))))
  cat("Proportion of non-default stepsize for xi =", 10, "is", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       observation = factor(rep(10, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     observation = factor(rep(10, nrepeats))))

load("inst/gaussian/results/integrator_y20.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       observation = factor(rep(20, total_non_default))))
  cat("Proportion of non-default stepsize for xi =", 20, "is", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       observation = factor(rep(20, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     observation = factor(rep(20, nrepeats))))

load("inst/gaussian/results/integrator_y30.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       observation = factor(rep(30, total_non_default))))
  cat("Proportion of non-default stepsize for xi =", 30, "is", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       observation = factor(rep(30, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     observation = factor(rep(30, nrepeats))))

load("inst/gaussian/results/integrator_y40.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       observation = factor(rep(40, total_non_default))))
  cat("Proportion of non-default stepsize for xi =", 40, "is", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       observation = factor(rep(40, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     observation = factor(rep(40, nrepeats))))

# plot no. of time steps against observation
gtimesteps <- ggplot(data = total.nsteps.df, aes(x = observation, y = timesteps)) + geom_boxplot() + 
  xlab(expression(xi)) + ylab("time steps")
gtimesteps
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_timesteps.eps", plot = gtimesteps,
       device = "eps", width = 6, height = 6)

# plot minimum step size against observation
gminstepsize <- ggplot(data = min.stepsize.df, aes(x = observation, y = minstepsize)) + geom_boxplot() + scale_y_log10() + 
  xlab(expression(xi)) + ylab("minimum stepsize")
gminstepsize
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_minstepsize.eps", plot = gminstepsize,
       device = "eps", width = 6, height = 6)
