\name{mvnplot}
\alias{mvnplot}
\title{Graphical and Numerical Tests of Multivariate Normality for Bootstrapped MDS Confidence Regions}
\description{This function implements the Friendly (1991) graphical test for multivariate normality on objects of class \code{bsmds}.  It also provides numerical summaries (p-values) for tests of multivariate normality based on skewness and kurtosis as suggested by Kankainen et. al. (2007)}
\usage{
mvnplot(obj, sim.ci=TRUE, nreps=1000, num=FALSE, plot=TRUE, method="integration")
}
\arguments{
  \item{obj}{An object of class \code{\link{bsmds}} that contains bootstrap replicates for a multidimensional scaling solution.}
  \item{sim.ci}{Logical indicating whether simulated confidence bounds should be placed around the qq-plot.}
  \item{nreps}{The number of simulations to preform to calculate the confidence bounds.}
  \item{num}{Logical indicating whether numerical tests of multinormality should be calculated and presented}
  \item{plot}{Logical indicating whether the plot is to be printed}
  \item{method}{Method for calculating p-values, must be one of \dQuote{integration}, \dQuote{satterthwaite} or dQuote{simulation}.  See \code{\link{mvnorm.kur.test}} from the \pkg{ICS} package for more details. }
}
\details{The Friendly (1991) graphical test plots the Mahalanobis distance between the bootstrap replications and the mean across the bootstrap replications against the quantiles of a chi-squared distribution with m (= number of dimensions in the multidimensional scaling solution) degrees of freedom.  If the bootstrap replicates follow a multivariate normal distribution, the points in the figure should fall on the 45-degree line which is imposed on the plot.  The details with respect to bootstrapped multidimensional scaling solutions are discussed in Jacoby and Armstrong (2011).  

The numerical tests presented are the p-values from the \code{\link{mvnorm.kur.test}} and \code{\link{mvnorm.skew.test}} functions from the \code{ICS} package as discussed in Kankainen et. a.. (2007).
}
\value{
\item{plot}{An \code{xyplot} object implementing the Friendly (1991) test}
\item{num}{A numeric matrix of p-values from the two numerical tests for multinormality discussed above.}
}
\references{
Friendly, M.  1991  \emph{SAS System for Statistical Graphics, 1st ed.}  SAS institute, Inc.

Jacoby, W.G. and D. Armstrong.  2011  Bootstrap Confidence Regions for Multidimensional Scaling Solutions.  Unpublished Manuscript.

Kankainen, A., S. Taskinen and H. Oja  2007  Tests of Multinormality Based on Location Vectors and Scatter Matrices, \emph{Statistical Methods and Applications} \bold{16}, 357--379.
}
\author{Dave Armstrong and William Jacoby}
\examples{
data(thermometers2004)
out <- bsmds(thermometers2004, dist.fun = "dist", dist.data.arg = "x", 
	dist.args=list(method="euclidian"), R=25)
mvnplot(out, sim.ci=TRUE, nreps=1000, num=TRUE)
}
\keyword{Statistics}
\keyword{Multivariate}





