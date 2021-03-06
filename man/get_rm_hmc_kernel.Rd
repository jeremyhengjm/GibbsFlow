% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rm_hmc_kernel.R
\name{get_rm_hmc_kernel}
\alias{get_rm_hmc_kernel}
\title{Construct Riemann manifold Hamiltonian Monte Carlo kernel for Cox process application}
\usage{
get_rm_hmc_kernel(logtarget, gradlogtarget, parameters)
}
\arguments{
\item{logtarget}{evaluates log density of the target distribution}

\item{gradlogtarget}{evaluates gradient of the log target density}

\item{parameters}{list with keys: 
\code{stepsize} specifies stepsize of leap-frog discretization, 
\code{nsteps} specifies number of leap-frog steps
\code{metric} specifies metric tensor}
}
\value{
list with key: \code{kernel}
}
\description{
Construct Riemann manifold Hamiltonian Monte Carlo kernel for Cox process application
}
