rm(list = ls())
library(ggplot2)
library(ggthemes)
nsteps <- 100
nrepeats <- 100
setmytheme()

# AIS dimension 4
load("inst/gaussian/results/repeat_ais_dim4.RData")
terminal_ess <- ess[, nsteps]
ess.df <- data.frame(lower = quantile(terminal_ess, probs = 0.25),
                     median = median(terminal_ess),
                     upper = quantile(terminal_ess, probs = 0.75),
                     method = "AIS", 
                     dimension = factor(4))

var_ais_dim4 <- var(log_normconst)
log_normconst.df <- data.frame(variance = var_ais_dim4, 
                               method = "AIS", 
                               dimension = factor(4))

# AIS dimension 8 
load("inst/gaussian/results/repeat_ais_dim8.RData")
terminal_ess <- ess[, nsteps]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "AIS", 
                                   dimension = factor(8)))

var_ais_dim8 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_ais_dim8, 
                                                       method = "AIS", 
                                                       dimension = factor(8)))

# Gibbs flow SIS dimension 4
load("inst/gaussian/results/repeat_gibbsflow_sis_dim4.RData")
terminal_ess <- ess[, nsteps]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "GF-SIS", 
                                   dimension = factor(4)))

var_gibbsflow_sis_dim4 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_gibbsflow_sis_dim4, 
                                                       method = "GF-SIS", 
                                                       dimension = factor(4)))
# Gibbs flow SIS dimension 8
load("inst/gaussian/results/repeat_gibbsflow_sis_dim8.RData")
terminal_ess <- ess[, nsteps]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "GF-SIS", 
                                   dimension = factor(8)))

var_gibbsflow_sis_dim8 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_gibbsflow_sis_dim8, 
                                                       method = "GF-SIS", 
                                                       dimension = 8))

# Gibbs flow AIS dimension 4
load("inst/gaussian/results/repeat_gibbsflow_ais_dim4.RData")
terminal_ess <- ess[, nsteps]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "GF-AIS", 
                                   dimension = factor(4)))

var_gibbsflow_ais_dim4 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_gibbsflow_ais_dim4, 
                                                       method = "GF-AIS", 
                                                       dimension = factor(4)))

# Gibbs flow AIS dimension 8
load("inst/gaussian/results/repeat_gibbsflow_ais_dim8.RData")
terminal_ess <- ess[, nsteps]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "GF-AIS", 
                                   dimension = factor(8)))

var_gibbsflow_ais_dim8 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_gibbsflow_ais_dim8, 
                                                       method = "GF-AIS", 
                                                       dimension = factor(8)))

# SVGD dimension 4
load("inst/gaussian/results/repeat_svgd_dim4.RData")
niterations <- ncol(ess)
terminal_ess <- ess[, niterations]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "SVGD", 
                                   dimension = factor(4)))

var_svgd_dim4 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_svgd_dim4, 
                                                       method = "SVGD", 
                                                       dimension = factor(4)))

# SVGD dimension 8
load("inst/gaussian/results/repeat_svgd_dim8.RData")
niterations <- ncol(ess)
terminal_ess <- ess[, niterations]
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "SVGD", 
                                   dimension = factor(8)))

var_svgd_dim8 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_svgd_dim8, 
                                                       method = "SVGD", 
                                                       dimension = factor(8)))

# TransportMaps (Marzouk et al) dimension 4
terminal_ess <- read.table(file = "inst/gaussian/results/transportmaps_ess_dim4.txt", sep = ",")
terminal_ess <- terminal_ess$V1
log_normconst <- read.table(file = "inst/gaussian/results/transportmaps_lognormconst_dim4.txt", sep = ",")
log_normconst <- log_normconst$V1
ess.df <- rbind(ess.df, data.frame(lower = quantile(terminal_ess, probs = 0.25),
                                   median = median(terminal_ess),
                                   upper = quantile(terminal_ess, probs = 0.75),
                                   method = "TM", 
                                   dimension = factor(4)))
var_transportmaps_dim4 <- var(log_normconst)
log_normconst.df <- rbind(log_normconst.df, data.frame(variance = var_transportmaps_dim4, 
                                                       method = "TM", 
                                                       dimension = factor(4)))

# ess plot
gess <- ggplot(ess.df, aes(x = method, y = median, ymin = lower, ymax = upper, color = dimension))
gess <- gess + geom_pointrange() + ylim(0, 100) + xlab("method") + ylab("ESS%") + scale_color_colorblind()
gess
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_ess.eps", plot = gess,
       device = "eps", width = 8, height = 6)

# variance of log-marginal likelihood 
remove_transportmaps <- !(log_normconst.df$method == "TM")
log_normconst.df <- log_normconst.df[remove_transportmaps, ]
gestimates <- ggplot(log_normconst.df, aes(x = method, y = variance, color = dimension))
gestimates <- gestimates + geom_point(size = 2) + scale_y_log10() + 
  xlab("method") + ylab(expression(Var(log(hat(Z)[M])))) + scale_color_colorblind() 
gestimates
ggsave(filename = "~/Dropbox/GibbsFlow/draft_v3/gaussian_estimate.eps", plot = gestimates,
       device = "eps", width = 7, height = 6)

# compare variance in dimension 4
cat("Variance of AIS:", var_ais_dim4)
cat("Variance of GF-SIS:", var_gibbsflow_sis_dim4)
cat("Variance of GF-AIS:", var_gibbsflow_ais_dim4)
cat("Variance of SVGD:", var_svgd_dim4)
cat("Variance of TransportMap:", var_transportmaps_dim4)

# compare variance in dimension 8
cat("Variance of AIS:", var_ais_dim8)
cat("Variance of GF-SIS:", var_gibbsflow_sis_dim8)
cat("Variance of GF-AIS:", var_gibbsflow_ais_dim8)
cat("Variance of SVGD:", var_gibbsflow_ais_dim8)

cat("Variance of SVGD / GF-AIS in dimension 4:", var_svgd_dim4 / var_gibbsflow_ais_dim4)
cat("Variance of SVGD / GF-AIS in dimension 8:", var_svgd_dim8 / var_gibbsflow_ais_dim8)
cat("Variance of SVGD / AIS in dimension 4:", var_svgd_dim4 / var_ais_dim4)
cat("Variance of SVGD / AIS in dimension 8:", var_svgd_dim8 / var_ais_dim8)

cat("Variance of AIS / GF-AIS in dimension 8:", var_ais_dim8 / var_gibbsflow_ais_dim8)
cat("Variance of AIS / GF-AIS:", var_ais / var_gibbsflow_ais)
cat("Variance of AIS / GF-SIS:", var_ais / var_gibbsflow_sis)
cat("Variance of GF-SIS / GF-AIS:", var_gibbsflow_sis / var_gibbsflow_ais)

