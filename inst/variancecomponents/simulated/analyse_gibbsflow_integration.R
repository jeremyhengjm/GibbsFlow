rm(list=ls())
library(ggplot2)
library(ggthemes)
setmytheme()

nrepeats <- 2^10
default_stepsize <- 0.01
exclude_default_stepsize <- TRUE

load("inst/variancecomponents/simulated/results/gibbsflow_integration_27.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- data.frame(minstepsize = min_stepsize[non_default], 
                                dimension = factor(rep(27, total_non_default)))
  cat("Dimension:", 27, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- data.frame(minstepsize = min_stepsize, 
                                dimension = factor(rep(27, nrepeats)))
}
total.nsteps.df <- data.frame(timesteps = total_nsteps, 
                              dimension = factor(rep(27, nrepeats)))

load("inst/variancecomponents/simulated/results/gibbsflow_integration_52.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(52, total_non_default))))
  cat("Dimension:", 52, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(52, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(52, nrepeats))))

load("inst/variancecomponents/simulated/results/gibbsflow_integration_102.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(102, total_non_default))))
  cat("Dimension:", 102, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(102, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(102, nrepeats))))

load("inst/variancecomponents/simulated/results/gibbsflow_integration_202.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(202, total_non_default))))
  cat("Dimension:", 202, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(202, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(202, nrepeats))))

load("inst/variancecomponents/simulated/results/gibbsflow_integration_402.RData")
if (exclude_default_stepsize){
  non_default <- (min_stepsize != default_stepsize)
  total_non_default <- sum(non_default)
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize[non_default], 
                                                       dimension = factor(rep(402, total_non_default))))
  cat("Dimension:", 402, "Percentage of non-default stepsizes:", total_non_default / nrepeats * 100, "%\n")
} else {
  min.stepsize.df <- rbind(min.stepsize.df, data.frame(minstepsize = min_stepsize, 
                                                       dimension = factor(rep(402, nrepeats))))
}
total.nsteps.df <- rbind(total.nsteps.df, data.frame(timesteps = total_nsteps, 
                                                     dimension = factor(rep(402, nrepeats))))

# plot no. of time steps against observation
gtimesteps <- ggplot(data = total.nsteps.df, aes(x = dimension, y = timesteps)) + geom_boxplot() + 
  xlab("dimension") + ylab("time steps")
gtimesteps
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_simulated_timesteps.eps", plot = gtimesteps,
       device = "eps", width = 6, height = 6)

# plot minimum step size against observation
gminstepsize <- ggplot(data = min.stepsize.df, aes(x = dimension, y = minstepsize)) + geom_boxplot() + scale_y_log10() + 
  xlab("dimension") + ylab("minimum stepsize")
gminstepsize
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_simulated_minstepsize.eps", plot = gminstepsize,
       device = "eps", width = 6, height = 6)

