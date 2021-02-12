rm(list=ls())
library(ggplot2)
library(ggthemes)
library(gridExtra)
setmytheme()

# results
load("inst/mixturemodel/smc_reference.RData")
load("inst/mixturemodel/gibbsflow_trajectory.RData")

# plot time evolution
marginals <- c(1, 2)
nsteps <- 400
times <- seq(0, 1, length.out = nsteps)
timesteps <- c(10, 50, 100, 400)
timings <- times[timesteps]
evolution_plots <- list()
for (i in 1:length(timesteps)){
  istep <- timesteps[i]
  
  xparticles <- smc$xtrajectory[ , , istep]
  reference_xparticles <- smc_reference$xtrajectory[ , , istep]
  plot.df <- data.frame(x = xparticles[, marginals[1]], y = xparticles[, marginals[2]], 
                        kde_x = reference_xparticles[, marginals[1]], kde_y = reference_xparticles[, marginals[2]])
  t <- as.character(signif(timings[i], 3))
  evolution_plots[[i]] <- ggplot(data = plot.df) +
    geom_density_2d(aes(x = kde_x, y = kde_y), bins = 5) +
    geom_point(aes(x = x, y = y), alpha = 1.0) +
    xlim(-10, 10) + 
    ylim(-10, 10) + 
    xlab(expression('x'[1])) + ylab(expression('x'[2])) + 
    labs(title = bquote(t == .(t)))
  
}
evolution_plot <- arrangeGrob(evolution_plots[[1]], evolution_plots[[2]], 
                              evolution_plots[[3]], evolution_plots[[4]],
                              ncol = 2, nrow = 2)
g <- grid.arrange(evolution_plot)
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_evolution.eps", plot = g,
       device = "eps", width = 8, height = 8)


# plot all pairwise marginals (at terminal time)
all_marginals <- combn(1:4,2) # 2 x 6 matrix
xparticles <- smc$xtrajectory[ , , 400]
reference_xparticles <- smc_reference$xtrajectory[ , , 400]
marginal_plots <- list()
for (k in 1:6){
  marginals <- all_marginals[, k]
  plot.df <- data.frame(x = xparticles[, marginals[1]], y = xparticles[, marginals[2]], 
                        kde_x = reference_xparticles[, marginals[1]], kde_y = reference_xparticles[, marginals[2]])
  index1 <- as.character(marginals[1])
  index2 <- as.character(marginals[2])
  marginal_plots[[k]] <- ggplot(data = plot.df) +
    geom_density_2d(aes(x = kde_x, y = kde_y), bins = 5) +
    geom_point(aes(x = x, y = y), alpha = 1.0) +
    xlim(-5, 8) + 
    ylim(-5, 8) + 
    xlab(bquote(x[.(index1)])) + 
    ylab(bquote(x[.(index2)]))  

}

marginal_plot <- arrangeGrob(marginal_plots[[1]], marginal_plots[[2]], 
                             marginal_plots[[3]], marginal_plots[[4]],
                             marginal_plots[[5]], marginal_plots[[6]],
                             ncol = 2, nrow = 3)
g <- grid.arrange(marginal_plot)
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_allpairs.eps", plot = g,
       device = "eps", width = 8, height = 12)

rm(list=ls())
library(ggthemes)
setmytheme()
load("inst/mixturemodel/gibbsflow_16384.RData")
load("inst/mixturemodel/plot_kde.RData")

marginals <- c(1, 2)
# plot kernel density estimates
plot.df <- data.frame(x = kde_smc$xparticles[, marginals[1]],
                      y = kde_smc$xparticles[, marginals[2]],
                      method = rep("SMC", 2^14))

plot.df <- rbind(plot.df, data.frame(x = smc$xparticles[, marginals[1]],
                                     y = smc$xparticles[, marginals[2]],
                                     method = rep("Gibbsflow", 2^14)))

g <- ggplot(data = plot.df) +
  geom_density_2d(aes(x = x, y = y, colour = method), bins = 5) +
  xlim(-5, 8) +
  ylim(-5, 8) + 
  xlab(expression('x'[1])) + ylab(expression('x'[2])) + 
  scale_colour_colorblind()
g 
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_gibbsflow_kdes.eps", plot = g,
       device = "eps", width = 6, height = 6)

# plot proportion of particles in each mode
true_means <- c(-3, 0, 3, 6)
nparticles <- nrow(smc$xparticles)
permutations <- function( x, prefix = c() )
{
  if(length(x) == 0 ) return(prefix)
  do.call(rbind, sapply(1:length(x), FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
}
cluster_particles <- kmeans(smc$xparticles, permutations(true_means))
cat("Empirical proportion:", table(cluster_particles$cluster) / nparticles)
cluster_particles$centers # check centers
barplot.df <- data.frame(index = cluster_particles$cluster)
g <- ggplot(barplot.df, aes(x = as.factor(index))) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  xlab("Mode") + 
  ylab("Proportion of samples")
g
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_gibbsflow_proportions.eps", plot = g,
       device = "eps", width = 10, height = 6)

# Pearson's Chi-squared test for uniformity 
test_output <- chisq.test(cluster_particles$size)
test_output$p.value

# results
rm(list=ls())
setmytheme()
load("inst/mixturemodel/gibbsflow_with_hmc.RData")
load("inst/mixturemodel/plot_kde.RData")

marginals <- c(1, 2)
# plot kernel density estimates
plot.df <- data.frame(x = kde_smc$xparticles[, marginals[1]],
                      y = kde_smc$xparticles[, marginals[2]],
                      method = rep("SMC", 2^14))

plot.df <- rbind(plot.df, data.frame(x = smc$xparticles[, marginals[1]],
                                     y = smc$xparticles[, marginals[2]],
                                     method = rep("Gibbsflow with HMC", 2^14)))

g <- ggplot(data = plot.df) +
  geom_density_2d(aes(x = x, y = y, colour = method), bins = 5) +
  xlim(-5, 8) +
  ylim(-5, 8) + 
  xlab(expression('x'[1])) + ylab(expression('x'[2])) + 
  scale_colour_colorblind()
g 
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_gibbsflow_hmc_kdes.eps", plot = g,
       device = "eps", width = 6, height = 6)
