rm(list = ls())
library(GibbsFlow)
library(ggplot2)
library(ggthemes)
setmytheme()

# posterior
posterior <- list()
posterior$cov <- function(lambda) gaussian_posterior_cov(lambda)

# (approximate) gibbs velocity field 
exponent <- 2
gibbs_velocity <- function(t, x){
  output <- gaussian_gibbsvelocity(t, x, exponent)
  return(output)
}

# true velocity field in 1D
true_velocity <- function(t, x){
  lambda <- t^exponent
  derivative_lambda <- exponent * t^(exponent - 1)
  output <- - derivative_lambda * x / (2 * (1 + lambda)) 
  return(output)
}

# compare approximate velocity field with true at t = 0.5
t <- 0.5
lambda <- t^exponent
current_std <- as.numeric(sqrt(posterior$cov(lambda)))
ngrid <- 100
range <- seq(-9 * current_std, 9 * current_std, length.out = ngrid)
std_breaks <- seq(-9, 9, by = 2) * current_std
std_labels <- c(expression(-9 * sigma), 
                expression(-7 * sigma),
                expression(-5 * sigma),
                expression(-3 * sigma),
                expression(-1 * sigma),
                expression(1 * sigma),
                expression(3 * sigma),
                expression(5 * sigma),
                expression(7 * sigma),
                expression(9 * sigma))
velocity.df <- data.frame(x = range, 
                          velocity = sapply(range, function(x) gibbs_velocity(t, x)),
                          type = factor(rep("Approximation", ngrid)))
velocity.df <- rbind(velocity.df, 
                     data.frame(x = range, 
                                velocity = true_velocity(t, range),
                                type = factor(rep("True", ngrid))))
g <- ggplot(data = velocity.df, aes(x = x, y = velocity, col = type)) + geom_line() + 
  scale_color_colorblind() + scale_x_continuous(breaks = std_breaks, 
                                                labels = std_labels) + 
  labs(x = expression(x))
g
pnorm(-8) 
pnorm(-9) 
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/plotvelocity1d_smallepsilon.eps", plot = g,
       device = "eps", width = 7, height = 6)
