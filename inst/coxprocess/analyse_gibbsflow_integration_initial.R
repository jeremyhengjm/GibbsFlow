rm(list=ls())
library(ggplot2)
library(ggthemes)
setmytheme()

nrepeats <- 2^10
default_stepsize <- 0.01
exclude_default_stepsize <- TRUE

# initialize from prior
load("inst/coxprocess/results/gibbsflow_integration_400.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- data.frame(minstepsize = min_stepsize[non_default], 
                                initialization = factor(rep("Prior", total_non_default)))
  cat("Initialization from prior:", "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- data.frame(minstepsize = min_stepsize, 
                                initialization = factor(rep("Prior", nrepeats)))
}
total.nsteps.df <- data.frame(timesteps = total_nsteps, 
                              initialization = factor(rep("Prior", nrepeats)))

# initialize from VB approximation 
load("inst/coxprocess/results/gibbsflow_integration_vi_400.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       initialization = factor(rep("VB", total_non_default))))
  cat("Initialization from VB:", "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       initialization = factor(rep("VB", nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     initialization = factor(rep("VB", nrepeats))))

# initialize from EP approximation 
load("inst/coxprocess/results/gibbsflow_integration_ep_400.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       initialization = factor(rep("EP", total_non_default))))
  cat("Initialization from EP:", "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       initialization = factor(rep("EP", nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     initialization = factor(rep("EP", nrepeats))))

# plot minimum step size against observation
gminstepsize <- ggplot(data = min.stepsize.df, aes(x = initialization, y = minstepsize)) + geom_boxplot() + scale_y_log10() + 
  xlab("initialization") + ylab("minimum stepsize")
gminstepsize
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_minstepsize_initial.eps", plot = gminstepsize,
       device = "eps", width = 6, height = 6)

# plot no. of time steps against observation
gtimesteps <- ggplot(data = total.nsteps.df, aes(x = initialization, y = timesteps)) + geom_boxplot() + 
  xlab("initialization") + ylab("time steps")
gtimesteps
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_timesteps_initial.eps", plot = gtimesteps,
       device = "eps", width = 6, height = 6)

