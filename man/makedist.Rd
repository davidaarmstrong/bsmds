\name{makedist}
\alias{makedist}
\title{Auxiliary function to make Dissimilarity Matrices}
\description{This function makes a distance matrix from a set of arguments}
\usage{
makedist(data, dist.fun, dist.data.arg, dist.args)
}
\arguments{
  \item{data}{A numeric data frame or matrix of variables to be scaled.  Column names will be propagated through the function as the names of stimuli.  This data frame must contain only the variables to be scaled.}
  \item{dist.fun}{A character string identifying the name of the distance function to be used.  We have included in this package a function to calculate Rabinowitz' Line-of-Sight distance.  The \code{dist} function from the MASS library is another option.  Users can provide their own function so long as it returns the distance/dissimilarity matrix in its object `dist'}
  \item{dist.data.arg}{A character string giving the name of the data argument to the distance function}
  \item{dist.args}{Optional arguments to be passed to the distance function - must not be NULL}
}
\value{
\item{dist}{Distance matrix returned from the distance function.  If the distance function is user-defined, it must return the distance matrix in its \code{dist} element.}
%\item{conv}{If the \code{los.dist} function is used to calculate distances, the \code{conv} element of the returned object indicates whether or not the algorithm converged.}
}
\author{Dave Armstrong and William Jacoby}
\examples{
data(thermometers2004)
d <- makedist(thermometers2004, dist.fun = "dist", dist.data.arg = "x", 
	dist.args=list(method="euclidian"))
d
}
\keyword{Statistics}
\keyword{Multivariate}
