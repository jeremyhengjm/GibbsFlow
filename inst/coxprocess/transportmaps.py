# import some stuff
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as npla
import scipy.stats as stats
import scipy.linalg as scila
import SpectralToolbox.Spectral1D as S1D
import TransportMaps as TM
import TransportMaps.Maps as MAPS
import TransportMaps.Distributions as DIST
import TransportMaps.Diagnostics as DIAG
import time
%matplotlib inline

# problem settings 
ngrid = 3;
grid = np.linspace(0, 1, ngrid+1)
dimension = ngrid * ngrid
parameter_sigmasq = 1.91
parameter_mu = np.log(126.0) - 0.5 * parameter_sigmasq
parameter_beta = 1.0 / 33.0
parameter_area = 1.0 / float(dimension)

# load pine saplings dataset
dataset = np.loadtxt("finpines.txt")
data_x = dataset[:,0]
data_y = dataset[:,1]
data_counts = np.zeros(dimension)
for i in range(ngrid):
    for j in range(ngrid):
        logical_y = np.logical_and(data_x > grid[i], data_x < grid[i+1])
        logical_x = np.logical_and(data_y > grid[j], data_y < grid[j+1])
        data_counts[i * ngrid + j] = float(sum(np.logical_and(logical_y, logical_x)))
        
# prior distribution
prior_mean = parameter_mu * np.ones(dimension)
prior_cov = np.zeros((dimension, dimension))
for m in range(dimension):
    for n in range(dimension):
        index_m = np.array([np.floor(m / ngrid) + 1, np.remainder(m, ngrid) + 1])
        index_n = np.array([np.floor(n / ngrid) + 1, np.remainder(n, ngrid) + 1])
        prior_cov[m, n] = parameter_sigmasq * np.exp(- np.sqrt(sum((index_m - index_n)**2)) / (float(ngrid) * parameter_beta))
prior = DIST.GaussianDistribution(mu=prior_mean, sigma=prior_cov)

# posterior distribution
class coxprocess(DIST.Distribution):
    def __init__(self):
        super(coxprocess, self).__init__(dimension)
    def pdf(self, x, params=None):
        return np.exp(self.log_pdf(x, params))
    #def log_pdf(self, x, params=None):
    def log_pdf(self, x, *args, **kwargs):
        out = prior.log_pdf(x) + np.sum(x * data_counts - parameter_area * np.exp(x), 1)
        return out.flatten()
    #def grad_x_log_pdf(self, x, params=None):
    def grad_x_log_pdf(self, x, *args, **kwargs):
        grad = prior.grad_x_log_pdf(x) + ( data_counts - parameter_area * np.exp(x) )
        return grad
    #def hess_x_log_pdf(self, x, params=None):
    def hess_x_log_pdf(self, x, *args, **kwargs):
        nrows = np.size(x, 0)
        hess = prior.hess_x_log_pdf(x)
        for n in range(nrows):
            hess[n,:,:] = - np.diag(x[n,:])
        return hess   
    
posterior = coxprocess()

# transport map settings
order = 1
T = TM.Default_IsotropicIntegratedSquaredTriangularTransportMap(dimension, order, 'full')
push_prior = DIST.PushForwardTransportMapDistribution(T, prior)
qtype = 3           # Gauss quadrature
qparams = [10] * dimension  # Quadrature order
reg = None          # No regularization
tol = 1e-5         # Optimization tolerance
ders = 2            # Use gradient and Hessian

# compute transport map
tStart = time.time()
optim_output = push_prior.minimize_kl_divergence(posterior, qtype=qtype, qparams=qparams, regularization=reg, tol=tol, ders=ders)
tEnd = time.time()
timeElapsed = tEnd - tStart
print("Time elapsed:", timeElapsed)
optim_output
