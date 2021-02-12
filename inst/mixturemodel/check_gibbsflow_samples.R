rm(list=ls())
library(ggplot2)
library(ggthemes)
library(gridExtra)
setmytheme()

# load reference
load("inst/mixturemodel/plot_kde.RData")
reference_xparticles <- kde_smc$xparticles

# load Gibbs flow samples
load("inst/mixturemodel/gibbsflow_16384.RData")
xparticles <- smc$xparticles

marginal_plots <- list()
marginals <- c(1, 2)
# plot Gibbs flow samples
plot.gibbs.df <- data.frame(x = xparticles[, marginals[1]], y = xparticles[, marginals[2]], 
                      kde_x = reference_xparticles[, marginals[1]], kde_y = reference_xparticles[, marginals[2]])
g <- ggplot(data = plot.gibbs.df) +
  geom_density_2d(aes(x = kde_x, y = kde_y), bins = 5) + # no. of contours = bins - 1
  geom_point(aes(x = x, y = y), alpha = 0.2, size = 2) +
  xlim(-5, 8) + 
  ylim(-5, 8) + 
  xlab(expression('x'[1])) + ylab(expression('x'[2])) + ggtitle('Gibbs flow')
g
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/mixturemodel_gibbsflow??.eps", plot = g,
       device = "eps", width = 6, height = 6)

# plot Gibbs flow samples
plot.smc.df <- data.frame(x = reference_xparticles[, marginals[1]], y = reference_xparticles[, marginals[2]], 
                          kde_x = reference_xparticles[, marginals[1]], kde_y = reference_xparticles[, marginals[2]])
g <- ggplot(data = plot.smc.df) +
  geom_density_2d(aes(x = kde_x, y = kde_y), bins = 5) +
  geom_point(aes(x = x, y = y), alpha = 0.2, size = 2) +
  xlim(-5, 8) + 
  ylim(-5, 8) + 
  xlab(expression('x'[1])) + ylab(expression('x'[2])) + ggtitle('Reference SMC')
g
# ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/mixturemodel_gibbsflow??.eps", plot = g,
# device = "eps", width = 6, height = 6)




