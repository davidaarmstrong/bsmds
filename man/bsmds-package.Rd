\name{bsmds-package}
\alias{bsmds-package}
\docType{package}
\title{Bootstrapping of Multidimensional Scaling Solutions}
\description{The \code{bsmds} package implements a bootstrapping solution for 
multidimensional scaling solutions based on the SMACOF function in the package
of the same name.}
\details{
\tabular{ll}{
Package: \tab bsmds\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2014-06-22\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The main function of the package is \code{bsmds} which implements the solution and returns the original configuration, the bootstrapped object (output from \code{boot}) and, for convenience, a number of matrices for diagnostic purposes. 
}
\author{
Dave Armstrong  and William Jacoby
Maintainer: Dave Armstrong <davearmstrong.ps@gmail.com>
}
\references{
Jacoby, William G. and David A. Armstrong II.  ``Bootstrap Confidence Regions for Multidimensional Scaling Solutions.''
}
\keyword{Statistics|Multivariate}
