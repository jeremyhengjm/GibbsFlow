rm(list=ls())
library(ggplot2)
library(ggthemes)
library(GibbsFlow)

M <- 80
nrepeats <- 100
setmytheme()
ess.df <- data.frame()
log_normconst.df <- data.frame() 

#### Prior initialization
# load AIS results
load("inst/coxprocess/results/repeat_ais_rmhmc_400.RData")
  
# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "Prior", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "AIS"))
  
# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "Prior",
                                                       variance = var(log_normconst), 
                                                       method = "AIS"))
var_ais <- sd(log_normconst)^2

# load GF-SIS results
load("inst/coxprocess/results/repeat_gibbsflow_sis_400.RData")
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "Prior", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-SIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "Prior",
                                                       variance = var(log_normconst), 
                                                       method = "GF-SIS"))
var_gibbsflow_sis <- sd(log_normconst)^2


# load GF-AIS results
load("inst/coxprocess/results/repeat_gibbsflow_ais_rmhmc_400.RData")
  
# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "Prior", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-AIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "Prior",
                                                       variance = var(log_normconst), 
                                                       method = "GF-AIS"))
var_gibbsflow_ais <- sd(log_normconst)^2

cat("Initialization: Prior", "\n")
cat("Variance of AIS:", var_ais, "\n")
cat("Variance of GF-SIS:", var_gibbsflow_sis, "\n")
cat("Variance of GF-AIS:", var_gibbsflow_ais, "\n")

cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais, "\n")
cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis, "\n")
cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais, "\n")

#### VB initialization
# load AIS results
load("inst/coxprocess/results/repeat_ais_vi_400.RData")

# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "VB", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "AIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "VB",
                                                       variance = var(log_normconst), 
                                                       method = "AIS"))
var_ais <- sd(log_normconst)^2

# load GF-SIS results
load("inst/coxprocess/results/repeat_gibbsflow_sis_vi_400.RData")
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "VB", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-SIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "VB",
                                                       variance = var(log_normconst), 
                                                       method = "GF-SIS"))
var_gibbsflow_sis <- sd(log_normconst)^2


# load GF-AIS results
load("inst/coxprocess/results/repeat_gibbsflow_ais_vi_400.RData")

# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "VB", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-AIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "VB",
                                                       variance = var(log_normconst), 
                                                       method = "GF-AIS"))
var_gibbsflow_ais <- sd(log_normconst)^2

cat("Initialization: VB", "\n")
cat("Variance of AIS:", var_ais, "\n")
cat("Variance of GF-SIS:", var_gibbsflow_sis, "\n")
cat("Variance of GF-AIS:", var_gibbsflow_ais, "\n")

cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais, "\n")
cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis, "\n")
cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais, "\n")

#### EP initialization
# load AIS results
load("inst/coxprocess/results/repeat_ais_ep_400.RData")

# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "EP", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "AIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "EP",
                                                       variance = var(log_normconst), 
                                                       method = "AIS"))
var_ais <- sd(log_normconst)^2

# load GF-SIS results
load("inst/coxprocess/results/repeat_gibbsflow_sis_ep_400.RData")
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "EP", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-SIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "EP",
                                                       variance = var(log_normconst), 
                                                       method = "GF-SIS"))
var_gibbsflow_sis <- sd(log_normconst)^2


# load GF-AIS results
load("inst/coxprocess/results/repeat_gibbsflow_ais_ep_400.RData")

# ess 
terminal_ess <- ess[, M]
ess.df <- rbind(ess.df, data.frame(initial = "EP", 
                                   lower = quantile(terminal_ess, probs = 0.25), 
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75), 
                                   method = "GF-AIS"))

# log-marginal likelihood
log_normconst.df <- rbind(log_normconst.df, data.frame(initial = "EP",
                                                       variance = var(log_normconst), 
                                                       method = "GF-AIS"))
var_gibbsflow_ais <- sd(log_normconst)^2

cat("Initialization: EP", "\n")
cat("Variance of AIS:", var_ais, "\n")
cat("Variance of GF-SIS:", var_gibbsflow_sis, "\n")
cat("Variance of GF-AIS:", var_gibbsflow_ais, "\n")

cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais, "\n")
cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis, "\n")
cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais, "\n")

# ess plot
gess <- ggplot(ess.df, aes(x = initial, y = median, ymin = lower, ymax = upper, color = method))
gess <- gess + geom_pointrange() + ylim(0, 100) + xlab("initialization") + ylab("ESS%") + scale_color_colorblind()
gess
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/coxprocess_initial_ess.eps", plot = gess,
       device = "eps", width = 6, height = 6)

# variance of log-marginal likelihood 
gestimates <- ggplot(log_normconst.df, aes(x = initial, y = variance, color = method))
gestimates <- gestimates + geom_point(size = 2) + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1.0)) + xlab("initialization") + ylab(expression(Var(log(hat(Z)[M])))) + scale_color_colorblind() 
gestimates
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v2/coxprocess_initial_estimate.eps", plot = gestimates,
       device = "eps", width = 6, height = 6)
