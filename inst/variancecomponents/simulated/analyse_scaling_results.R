rm(list=ls())
library(ggplot2)
library(ggthemes)
library(GibbsFlow)
setmytheme()

dimensions <- c(27, 52, 102, 202, 402)
ndimensions <- length(dimensions)
nsteps <- c(125, 250, 500, 750, 1000)
nrepeats <- 100
setmytheme()
ess.df <- data.frame()
log_normconst.df <- data.frame() 

for (i in 1:ndimensions){
  d <- dimensions[i]
  M <- nsteps[i]
  
  # load AIS results
  filename = paste("inst/variancecomponents/simulated/results/repeat_ais_", d, ".RData", sep = "")
  load(file = filename)
  
  # ess 
  terminal_ess <- ess[, M]
  ess.df <- rbind(ess.df, data.frame(dimension = d, 
                                     lower = quantile(terminal_ess, probs = 0.25), 
                                     median = median(terminal_ess),
                                     upper = quantile(terminal_ess, probs = 0.75), 
                                     method = "AIS"))
  
  # log-marginal likelihood
  log_normconst.df <- rbind(log_normconst.df, data.frame(dimension = d,
                                                         variance = var(log_normconst), 
                                                         method = "AIS"))
  var_ais <- sd(log_normconst)^2
  
  # load GF-SIS results
  filename = paste("inst/variancecomponents/simulated/results/repeat_gibbsflow_sis_", d, ".RData", sep = "")
  load(file = filename)
  
  # ess 
  terminal_ess <- ess[, M]
  ess.df <- rbind(ess.df, data.frame(dimension = d, 
                                     lower = quantile(terminal_ess, probs = 0.25), 
                                     median = median(terminal_ess),
                                     upper = quantile(terminal_ess, probs = 0.75), 
                                     method = "GF-SIS"))
  
  # log-marginal likelihood
  log_normconst.df <- rbind(log_normconst.df, data.frame(dimension = d,
                                                         variance = var(log_normconst), 
                                                         method = "GF-SIS"))
  var_gibbsflow_sis <- sd(log_normconst)^2
  
  # load GF-AIS results
  filename = paste("inst/variancecomponents/simulated/results/repeat_gibbsflow_ais_", d, ".RData", sep = "")
  load(file = filename)
  
  # ess 
  terminal_ess <- ess[, M]
  ess.df <- rbind(ess.df, data.frame(dimension = d, 
                                     lower = quantile(terminal_ess, probs = 0.25), 
                                     median = median(terminal_ess),
                                     upper = quantile(terminal_ess, probs = 0.75), 
                                     method = "GF-AIS"))
  
  # log-marginal likelihood
  log_normconst.df <- rbind(log_normconst.df, data.frame(dimension = d,
                                                         variance = var(log_normconst), 
                                                         method = "GF-AIS"))
  var_gibbsflow_ais <- sd(log_normconst)^2
  
  cat("Dimension:", d, "\n")
  cat("Variance of AIS:", var_ais, "\n")
  cat("Variance of GF-SIS:", var_gibbsflow_sis, "\n")
  cat("Variance of GF-AIS:", var_gibbsflow_ais, "\n")
  
  cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais, "\n")
  cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis, "\n")
  cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais, "\n")
  
  
}

# ess plot
gess <- ggplot(ess.df, aes(x = dimension, y = median, ymin = lower, ymax = upper, color = method))
gess <- gess + geom_pointrange() + xlim(1,500) + ylim(0, 100) + xlab("dimension") + ylab("ESS%") + scale_color_colorblind()
gess
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/vcmodel_simulated_ess.eps", plot = gess,
       device = "eps", width = 6, height = 6)

# variance of log-marginal likelihood 
gestimates <- ggplot(log_normconst.df, aes(x = dimension, y = variance, color = method))
gestimates <- gestimates + geom_point(size = 2) + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1.0)) + xlim(1,500) + 
              xlab("dimension") + ylab(expression(Var(log(hat(Z)[M])))) + scale_color_colorblind() 
gestimates
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/vcmodel_simulated_estimate.eps", plot = gestimates,
       device = "eps", width = 6, height = 6)


