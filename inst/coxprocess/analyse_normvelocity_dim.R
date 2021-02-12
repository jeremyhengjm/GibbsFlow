library(ggplot2)
library(ggthemes)
library(GibbsFlow)
setmytheme()

# plot norm of Gibbs velocity for GF-SIS for different dimensions
rm(list=ls())
nsteps <- 40
load("inst/coxprocess/results/normvelocity_gfsis_100.RData")
plot.normvelocity.df <- cbind(normvelocity.df, data.frame(dimension = factor(rep("100", nsteps))))

nsteps <- 60
load("inst/coxprocess/results/normvelocity_gfsis_225.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = factor(rep("225", nsteps)))))
                              
nsteps <- 80
load("inst/coxprocess/results/normvelocity_gfsis_400.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("400", nsteps))))

gnormvelocity <- ggplot(plot.normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper, colour = dimension))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 1) + 
  xlim(0, 1) + ylim(0, 20) + 
  # scale_y_continuous(breaks = c(0, 40, 80, 120)) +
  xlab("time") + ylab("norm of Gibbs velocity") + scale_color_colorblind()
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_normvelocity_gfsis_dim.eps", plot = gnormvelocity,
       device = "eps", width = 6, height = 6)


# plot norm of Gibbs velocity for GF-AIS for different dimensions
rm(list=ls())
nsteps <- 40
load("inst/coxprocess/results/normvelocity_gfais_100.RData")
plot.normvelocity.df <- cbind(normvelocity.df, data.frame(dimension = factor(rep("100", nsteps))))

nsteps <- 60
load("inst/coxprocess/results/normvelocity_gfais_225.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = factor(rep("225", nsteps)))))

nsteps <- 80
load("inst/coxprocess/results/normvelocity_gfais_400.RData")
plot.normvelocity.df <- rbind(plot.normvelocity.df, 
                              cbind(normvelocity.df, data.frame(dimension = rep("400", nsteps))))

gnormvelocity <- ggplot(plot.normvelocity.df, aes(x = time, y = median, ymin = lower, ymax = upper, colour = dimension))
gnormvelocity <- gnormvelocity + geom_pointrange(alpha = 1) + 
  xlim(0, 1) + ylim(0, 20) + 
  # scale_y_continuous(breaks = c(0, 40, 80, 120)) +
  xlab("time") + ylab("norm of Gibbs velocity") + scale_color_colorblind()
gnormvelocity
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/coxprocess_normvelocity_gfais_dim.eps", plot = gnormvelocity,
device = "eps", width = 6, height = 6)