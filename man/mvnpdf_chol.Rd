% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mvnpdf_chol}
\alias{mvnpdf_chol}
\title{Evaluate multivariate Gaussian density (with pre-computed Cholesky factor)}
\usage{
mvnpdf_chol(x, mu, rooti, rootisum)
}
\arguments{
\item{x}{evaluation points}

\item{mu}{mean vector}

\item{rooti}{inverse of Cholesky factor}

\item{rootisum}{normalizing constant term}
}
\value{
density values
}
\description{
Evaluate multivariate Gaussian density (with pre-computed Cholesky factor)
}
