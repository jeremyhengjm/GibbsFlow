rm(list = ls())
library(ggplot2)
library(ggthemes)
nsteps <- 50
nrepeats <- 100
setmytheme()

load("inst/variancecomponents/baseball/results/repeat_ais.RData")

ess.df <- data.frame(time = seq(0, 1, length.out = nsteps), 
                     lower = apply(ess, 2, function(x) quantile(x, probs = 0.25)),
                     median = apply(ess, 2, median),
                     upper = apply(ess, 2, function(x) quantile(x, probs = 0.75)),
                     method = rep("AIS", nsteps))

log_normconst.df <- data.frame(estimates = log_normconst, method = rep("AIS", nrepeats))
var_ais <- sd(log_normconst)^2

load("inst/variancecomponents/baseball/results/repeat_gibbsflow_sis.RData")

ess.df <- rbind(ess.df, data.frame(time = seq(0, 1, length.out = nsteps), 
                                   lower = apply(ess, 2, function(x) quantile(x, probs = 0.25)),
                                   median = apply(ess, 2, median),
                                   upper = apply(ess, 2, function(x) quantile(x, probs = 0.75)),
                                   method = rep("GF-SIS", nsteps)))

log_normconst.df <- rbind(log_normconst.df, data.frame(estimates = log_normconst, 
                                                       method = rep("GF-SIS", nrepeats)))
var_gibbsflow_sis <- sd(log_normconst)^2

load("inst/variancecomponents/baseball/results/repeat_gibbsflow_ais.RData")

ess.df <- rbind(ess.df, data.frame(time = seq(0, 1, length.out = nsteps), 
                                   lower = apply(ess, 2, function(x) quantile(x, probs = 0.25)),
                                   median = apply(ess, 2, median),
                                   upper = apply(ess, 2, function(x) quantile(x, probs = 0.75)),
                                   method = rep("GF-AIS", nsteps)))

log_normconst.df <- rbind(log_normconst.df, data.frame(estimates = log_normconst, 
                                                       method = rep("GF-AIS", nrepeats)))
cat("Std of Gibbs flow AIS estimates:", sd(log_normconst))
var_gibbsflow_ais <- sd(log_normconst)^2
  
gess <- ggplot(ess.df, aes(x = time, y = median, ymin = lower, ymax = upper, color = method))
gess <- gess + geom_pointrange() + xlim(0, 1) + ylim(0, 100) + xlab("time") + ylab("ESS%") + scale_color_colorblind()
gess
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/vcmodel_baseball_ess.eps", plot = gess,
       device = "eps", width = 6, height = 6)

gestimates <- ggplot(log_normconst.df, aes(x = method, y = estimates))
gestimates <- gestimates + geom_boxplot() + ylab(expression(log(hat(Z)[M])))
gestimates
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/vcmodel_baseball_estimate.eps", plot = gestimates,
       device = "eps", width = 6, height = 6)

cat("Variance of AIS:", var_ais)
cat("Variance of GF-SIS:", var_gibbsflow_sis)
cat("Variance of GF-AIS:", var_gibbsflow_ais)

cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais)
cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis)
cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais)

