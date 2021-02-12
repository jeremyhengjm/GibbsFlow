rm(list=ls())
library(ggplot2)
library(ggthemes)
library(GibbsFlow)
setmytheme()

# plot norm of Gibbs velocity for GF-SIS for different dimensions
nsteps <- 125
load("inst/variancecomponents/simulated/results/normvelocity_gfsis_27.RData")
plot.normvelocity.df <- cbind(normvelocity.df, data.frame(dimension = factor(rep("27", nsteps))))

nsteps <- 250
load("inst/variancecomponents/simulated/results/normvelocity_gfsis_52.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = factor(rep("52", nsteps)))))

nsteps <- 500
load("inst/variancecomponents/simulated/results/normvelocity_gfsis_102.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("102", nsteps))))

nsteps <- 750
load("inst/variancecomponents/simulated/results/normvelocity_gfsis_202.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("202", nsteps))))

nsteps <- 1000
load("inst/variancecomponents/simulated/results/normvelocity_gfsis_402.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("402", nsteps))))

gnormvelocity <- ggplot(plot.normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper, colour = dimension))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 0.5) + 
  xlim(0, 1) + ylim(0, 10) + 
  xlab("time") + ylab("norm of Gibbs velocity") + scale_color_colorblind()
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_simulated_normvelocity_gfsis_dim.pdf", plot = gnormvelocity,
       device = "pdf", width = 8, height = 6)

# plot norm of Gibbs velocity for GF-AIS for different dimensions
rm(list=ls())
nsteps <- 125
load("inst/variancecomponents/simulated/results/normvelocity_gfais_27.RData")
plot.normvelocity.df <- cbind(normvelocity.df, data.frame(dimension = factor(rep("27", nsteps))))

nsteps <- 250
load("inst/variancecomponents/simulated/results/normvelocity_gfais_52.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = factor(rep("52", nsteps)))))

nsteps <- 500
load("inst/variancecomponents/simulated/results/normvelocity_gfais_102.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("102", nsteps))))

nsteps <- 750
load("inst/variancecomponents/simulated/results/normvelocity_gfais_202.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("202", nsteps))))

nsteps <- 1000
load("inst/variancecomponents/simulated/results/normvelocity_gfais_402.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("402", nsteps))))

gnormvelocity <- ggplot(plot.normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper, colour = dimension))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 0.5) + 
  xlim(0, 1) + ylim(0, 10) + 
  xlab("time") + ylab("norm of Gibbs velocity") + scale_color_colorblind()
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/vcmodel_simulated_normvelocity_gfais_dim.pdf", plot = gnormvelocity,
       device = "pdf", width = 8, height = 6)
