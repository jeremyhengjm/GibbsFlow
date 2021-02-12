rm(list=ls())
library(ggplot2)
library(ggthemes)
# library(gridExtra)
setmytheme()

plot_covariance <- function(matrix, npoints){
  # compute eigendecomposition
  eigens <- eigen(matrix)
  evs <- sqrt(eigens$values)
  evecs <- eigens$vectors
  
  a <- evs[1]
  b <- evs[2]
  x0 <- 0
  y0 <- 0
  alpha <- atan(evecs[ , 1][2] / evecs[ , 1][1])
  theta <- seq(0, 2 * pi, length = npoints)
  
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
  
  return(list(x = x, y = y))
}

# select marginals
marginals <- c(1, 2)

# load reference
load("inst/mixturemodel/plot_kde.RData")
reference_xparticles <- kde_smc$xparticles[, marginals]

# load Gibbs flow samples
load("inst/mixturemodel/gibbsflow_16384.RData")
gibbs_xparticles <- smc$xparticles[, marginals]

# cluster
nmodes <- 4 * 3
true_means <- c(-3, 0, 3, 6)
modes <- matrix(0, nrow = nmodes, ncol = 2)
modes[1:6, ] <- t(combn(true_means, 2))
modes[7:12, 1] <- modes[1:6, 2]
modes[7:12, 2] <- modes[1:6, 1]

cluster_reference <- kmeans(reference_xparticles, centers = modes)
cluster_gibbs <- kmeans(gibbs_xparticles, centers = modes)

# compute sample covariance matrix in each mode 
plot.df <- data.frame()
npoints <- 100
for (i in 1:nmodes){
  
  # gibbs particles
  members <- (cluster_gibbs$cluster == i)
  member_particles <- gibbs_xparticles[members, ]
  mean_vector <- colMeans(member_particles)
  covariance_matrix <- cov(member_particles)
  plot_matrix <- plot_covariance(covariance_matrix, npoints)
  plot_matrix$x <- plot_matrix$x + mean_vector[1]
  plot_matrix$y <- plot_matrix$y + mean_vector[2]
  plot.df <- rbind(plot.df, data.frame(x = plot_matrix$x,
                                       y = plot_matrix$y,
                                       method = rep("Gibbsflow", npoints)))
  
  # reference particles
  members <- (cluster_reference$cluster == i)
  member_particles <- reference_xparticles[members, ]
  mean_vector <- colMeans(member_particles)
  covariance_matrix <- cov(member_particles)
  plot_matrix <- plot_covariance(covariance_matrix, npoints)
  plot_matrix$x <- plot_matrix$x + mean_vector[1]
  plot_matrix$y <- plot_matrix$y + mean_vector[2]
  plot.df <- rbind(plot.df, data.frame(x = plot_matrix$x, 
                                       y = plot_matrix$y, 
                                       method = rep("SMC", npoints)))
                   
}

g <- ggplot(data = plot.df) + geom_point(aes(x = x , y = y, col = method), size = 0.1) + 
  xlim(-5, 8) +
  ylim(-5, 8) +
  xlab(expression('x'[1])) + ylab(expression('x'[2])) + 
  scale_colour_colorblind() 
g

ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_compare_cov_marginal12.eps", plot = g,
       device = "eps", width = 6, height = 6)

# select marginals
marginals <- c(3, 4)
reference_xparticles <- kde_smc$xparticles[, marginals]
gibbs_xparticles <- smc$xparticles[, marginals]

# cluster
nmodes <- 4 * 3
true_means <- c(-3, 0, 3, 6)
modes <- matrix(0, nrow = nmodes, ncol = 2)
modes[1:6, ] <- t(combn(true_means, 2))
modes[7:12, 1] <- modes[1:6, 2]
modes[7:12, 2] <- modes[1:6, 1]

cluster_reference <- kmeans(reference_xparticles, centers = modes)
cluster_gibbs <- kmeans(gibbs_xparticles, centers = modes)

# compute sample covariance matrix in each mode 
plot.df <- data.frame()
npoints <- 100
for (i in 1:nmodes){
                   
  # gibbs particles
  members <- (cluster_gibbs$cluster == i)
  member_particles <- gibbs_xparticles[members, ]
  mean_vector <- colMeans(member_particles)
  covariance_matrix <- cov(member_particles)
  plot_matrix <- plot_covariance(covariance_matrix, npoints)
  plot_matrix$x <- plot_matrix$x + mean_vector[1]
  plot_matrix$y <- plot_matrix$y + mean_vector[2]
  plot.df <- rbind(plot.df, data.frame(x = plot_matrix$x,
                                       y = plot_matrix$y,
                                       method = rep("Gibbsflow", npoints)))
  
  # reference particles
  members <- (cluster_reference$cluster == i)
  member_particles <- reference_xparticles[members, ]
  mean_vector <- colMeans(member_particles)
  covariance_matrix <- cov(member_particles)
  plot_matrix <- plot_covariance(covariance_matrix, npoints)
  plot_matrix$x <- plot_matrix$x + mean_vector[1]
  plot_matrix$y <- plot_matrix$y + mean_vector[2]
  plot.df <- rbind(plot.df, data.frame(x = plot_matrix$x, 
                                       y = plot_matrix$y, 
                                       method = rep("SMC", npoints)))
  
}

g <- ggplot(data = plot.df) + geom_point(aes(x = x , y = y, col = method), size = 0.1) + 
  xlim(-5, 8) +
  ylim(-5, 8) +
  xlab(expression('x'[3])) + ylab(expression('x'[4])) + 
  scale_colour_colorblind() 
g

ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/mixturemodel_compare_cov_marginal34.eps", plot = g,
       device = "eps", width = 6, height = 6)
