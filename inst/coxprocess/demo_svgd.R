rm(list = ls())
library(GibbsFlow)
library(spatstat)
library(tictoc)
library(ggplot2)

# load pine saplings dataset
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
plot(x = data_x, y = data_y, type = "p")

ngrid <- 20
grid <- seq(from = 0, to = 1, length.out = ngrid+1)
dimension <- ngrid^2
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}

# Stein variational importance sampling (approx 28s)
nfparticles <- 32
nlparticles <- 40
niterations <- 10
stepsize <- 0.1
tic()
  svis <- coxprocess_stein_variational_importance_sampling(nfparticles, nlparticles, niterations, stepsize, data_counts)
toc()

svis$ess
svis$log_normconst

# repeats
nfparticles <- 32
nlparticles <- 40
niterations <- 10
stepsize <- 0.1

nrepeats <- 100
ess <- rep(0, nrepeats)
log_normconst <- rep(0, nrepeats)

tic()
for (i in 1:nrepeats){
  cat("Repeat: ", i, " / ", nrepeats, "\n")
  svis <- coxprocess_stein_variational_importance_sampling(nfparticles, nlparticles, niterations, stepsize, data_counts)
  ess[i] <- svis$ess
  log_normconst[i] <- svis$log_normconst
}
toc()
save(ess, log_normconst, file = "inst/coxprocess/results/repeat_svis.RData", safe = F)

cat("Average ESS%", mean(ess), "\n")
cat("Variance log normalizing constant estimator:", var(log_normconst), "\n")

