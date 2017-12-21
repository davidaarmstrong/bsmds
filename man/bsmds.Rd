\name{bsmds}
\alias{bsmds}
\alias{bsfun}
\title{Bootstrap Confidence Regions for Multidimensional Scaling Solutions}
\description{This function implements a bootstrapping solution that identifies and optionally removes potentially degenerate solutions.}
\usage{
bsmds(data, dist.fun, dist.data.arg = "x", dist.args = NULL, R, ndim = 2, weightmat = NULL, init = NULL, type="interval", ties = "primary", verbose = FALSE, relax = 1, modulus = 1, itmax = 1000, eps = 1e-06, rm.degen = TRUE, km.thresh = 5, iter.info=FALSE)
}
\arguments{
  \item{data}{A numeric data frame or matrix of variables to be scaled.  Column names will be propagated through the function as the names of stimuli.  This data frame must contain only the variables to be scaled.}
  \item{dist.fun}{A character string identifying the name of the distance function to be used.  We have included in this package a function to calculate Rabinowitz' Line-of-Sight distance.  The \code{dist} function from the MASS library is another option.  Users can provide their own function so long as it returns the distance/dissimilarity matrix in its object `dist'}
  \item{dist.data.arg}{A character string giving the name of the data argument to the distance function}
  \item{dist.args}{Optional arguments to be passed to the distance function - must not be NULL}
  \item{R}{Number of bootstrap replications to compute}
  \item{ndim}{Number of dimensions for the multidimensional scaling solution}
  \item{weightmat}{Optional matrix of dissimilarity weights to be passed to \code{smacofSym}}
  \item{init}{Matrix with starting values for the initial configuration to be passed to \code{smacofSym}, in lieu of a matrix, supply \sQuote{togrgerson} (default) or \sQuote{random}.}
  \item{type}{How should distances be treated (formerly the \code{METRIC} argument), can be \sQuote{interval}, \sQuote{ratio}, \sQuote{mspline}, \sQuote{ordinal}.}
  \item{ties}{Tie specification for non-metric MDS only: \code{"primary"}, \code{"secondary"}, or \code{"tertiary"}, to be passed to \code{smacofSym}}
  \item{verbose}{If \code{TRUE}, intermediate stress is printed out}
  \item{relax}{If \code{TRUE}, block relaxation is used for majorization, to be passed to \code{smacofSym}}
  \item{modulus}{Number of smacof iterations per monotone regression call}
  \item{itmax}{Maximum number of iterations}
  \item{eps}{Convergence criterion}
  \item{rm.degen}{A logical argument indicating whether the algorithm should remove potentially degenerate solutions.  See `Details' for more information}
  \item{km.thresh}{Integer-valued threshold (between 2 and n-stimuli-1) of the maximum number of clusters for which, if approximately 100\% of the total variance is accounted for by between-cluster variance, the solution is deemed to be degenerate}
  \item{iter.info}{Logical indicating whether information about each iteration should be printed}
}
\details{The \code{bsmds} function bootstraps data to create bootstrap replicates of the dissimilarity matrix on which the MDS solution can be computed with \code{smacofSym}.   Our experience is that this procedure has a tendency to generate a non-trivial proportion of degenerate solutions (i.e., those where the scaled stimuli locations fall into a small number of clusters).  We include an optional k-means clustering step (which can be turned on by \code{rm.degen = TRUE}) which calculates the ratio of between-cluster to total sums of squares.  If this ratio is approximately 1 for any number of clusters less than or equal to \code{km.thresh}, the solution is deemed to be degenerate and a different bootstrap sample is drawn and a new dissimilarity matrix is calculated.\\  

The function performs a Procrustean similarity transformation, proposed by Schonemann and Carroll (1970), to the bootstrap configuration that brings each into maximum similarity with the target configuration (the configuration from the original MDS solution).  This consists of a translation, rotation and dilation of the bootstrap configuration.
 }
\value{The returned object is of class \code{bsmds} with elements:
\item{mds}{The original mds solution returned by \code{smacofSym}}
\item{bs.mds}{The bootstrapped solution (the return from a call to \code{boot})}
\item{dmat}{The original dissimilarities matrix used to produce the scaling solution}
\item{X.i}{A list of bootstrap stimuli locations.  Each element of the list is an R by n-dim matrix of points for each stimulus.  The names of the list elements come from the column names of the data matrix}
\item{cors}{An R by n-dim*2 matrix of correlations.  The first n-dim columns give the correlations between the original MDS configuration and the untransformed bootstrap configuration.  The last n-dim columns of the matrix are the correlations between the optimally transformed bootstrapped configurations and the original MDS configuration}
\item{v}{A vector of values indicating whether the solution was valid (1) or degenerate (0)}
\item{stress}{The stress (either metric or non-metric, as appropriate) from each bootstrap replication of the \code{smacofSym} procedure}
\item{pct.totss}{The percentage of the total sum of squares accounted for by the between-cluster sum of squares in the k-means clustering solution.  This is returned regardless of whether \code{rm.degen=TRUE}, though solutions are only discarded on this bases if that argument is turned on.}}
\references{
de Leeuw, J. & Mair, P. (2009). Multidimensional scaling using majorization: 
The R package smacof. Journal of Statistical Software, 31(3), 1-30, \url{http://www.jstatsoft.org/v31/i03/}

Jacoby, W.G. and D. Armstrong.  2011  Bootstrap Confidence Regions for Multidimensional Scaling Solutions.  Unpublished Manuscript.

Schonemann, P. and R. M. Carroll.  1970  \dQuote{Fitting One Matrix to Another Under Choice of a Central Dilation and a Rigid Motion} \emph{Psychometrika} \bold{35}, 245--256.
}
\author{Dave Armstrong and William Jacoby}
\examples{
data(thermometers2004)
out <- bsmds(thermometers2004, dist.fun = "dist", dist.data.arg = "x", 
	dist.args=list(method="euclidian"), R=25, type="interval", iter.info=FALSE)
}
\keyword{Statistics}
\keyword{Multivariate}
