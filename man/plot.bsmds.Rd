\name{plot.bsmds}
\alias{plot.bsmds}
\title{Plot method for objects of class bsmds}
\description{This method plots either the configuration with bootstrap confidence ellipses or a set of diagnostic measures to investigate the solution.}
\usage{
\method{plot}{bsmds}(x, ..., parametric = TRUE, col = NULL, lwd = 1, center.cex = 1, id = c("colors", "text", "number", "pch", "none"), confidence.level = 0.95, key.side = "top", key.ncol = 3, srt = 0, loc = NULL, interactive.label.placement = FALSE, text.cex = 1, type = "config")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{bsmds}}
  \item{...}{Optional arguments to be passed to \code{plot}.}
  \item{parametric}{Logical argument whether normal-theory confidence ellipses are plotted (if \code{TRUE}) or convex hulls are drawn (if \code{FALSE}).}
  \item{col}{A vector of colors used to identify the different ellipses.}
  \item{lwd}{Width of the lines used to draw the ellipses or hulls.}
  \item{center.cex}{The size of the point in the middle of the ellipse or hull.}
  \item{id}{String indicating how the ellipses should be identified.  There are a number of options including \code{"colors"} in which case the ellipses are drawn with an optionally user-provided vector of colors, \code{"text"} in which case the stimulus identifiers will be plotted in the middle of the ellipses, \code{"number"} in which case a number is printed in the middle of the ellipse and identified in the key to the plot, \code{"pch"} in which case different plotting symbols are used and \code{"none"} in which case no identifier is used.}
  \item{confidence.level}{The confidence level of the ellipse or hull to be plotted.}
  \item{key.side}{The side on which the key should to the graph should be placed.}
  \item{key.ncol}{The number of columns in the key.}
  \item{srt}{An optional vector indicating the number of degrees each label should be rotated.}
  \item{loc}{An optional vector of locations for the identifying text to be plotted.}
  \item{interactive.label.placement}{If set to \code{TRUE}, allows the user to interactively (re)position the text labels for the ellipses.}
  \item{text.cex}{An optional vector of values that indicates the character expansion for the text identifiers.}
  \item{type}{If \code{"config"}, the configuration with bootstrap confidence regions is printed.  If \code{"diag"} a series of diagnostics are plotted.  See \sQuote{Details} for more information.}
}
\details{
When parametric confidence regions are plotted, the \code{ellipse} function from the package of the same name is used.  Otherwise, when non-parametric regions are desired, the \code{chull} function from the \code{grDevices} package is used to find the convex hull for the 100(\code{confidence.level})\% of points with the smallest euclidian distance to the original configuration point. 

The interactive label placement feature provides a menu-driven interface where users can place the label and then optionally reposition it if they choose.  If \code{loc} is \code{NULL} and \code{interactive.label.placement} is \code{FALSE}, the \code{pointLabel} function from the \code{maptools} package is used to find optimal placement of the labels with the smallest amount of overlap.

If the user requests \code{type = "diag"}, then three diagnostic plots are printed.  First, the between-cluster SS as a proportion of total SS from 1:\code{km.thresh} are printed for each bootstrap sample.  Those iterations identified as degenerate are plotted in red.  Second, a box-plot of stress (either metric or non-metric, as appropriate) is plotted for both the valid and degenerate solutions.  Finally, a scatterplot of correlations is printed where the x-axis represents the correlation between the target (the original MDS configuration) and the un-transformed bootstrap configuration and the y-axis represents the correlation between the target and the optimally transformed bootstrap configuration.  The plotting symbols are the dimension numbers.  The values used to make the plot are in the \code{cors} element of the \code{bsmds} object if further inspection is required. 

}
\value{
\item{plot}{The \code{lattice} object is returned if \code{type = "config"}}
\item{loc}{If \code{type = "config"}, then the location of the text labels (if appropriate) is returned}
}
\references{
Jacoby, W.G. and D. Armstrong.  2011  Bootstrap Confidence Regions for Multidimensional Scaling Solutions.  Unpublished Manuscript.
}
\author{Dave Armstrong and William Jacoby}
\examples{
data(thermometers2004)
out <- bsmds(thermometers2004, dist.fun = "dist", dist.data.arg="x",
    dist.args=list(method="euclidian"), R=25)
plot(out, type="diag")
plot(out, type="config", id="number")
}
\keyword{Graphics}
\keyword{Statistics}
\keyword{Multivariate}
