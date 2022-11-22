\name{persp.polymars}
\alias{persp.polymars}
\title{Polymars: multivariate adaptive polynomial spline regression}
\description{This function is not intended for direct use. It is called
by \code{\link{plot.polymars}}.}
\usage{\method{persp}{polymars}(x, predictor1, predictor2, response, n = 33,
xlim, ylim, xx, contour.polymars, main, intercept, ...) }
\arguments{
\item{x, predictor1, predictor2}{this function is not intended to be called directly.}
\item{response, n, xlim, ylim}{this function is not intended to be called directly.}
\item{xx, contour.polymars}{this function is not intended to be called directly.}
\item{main, intercept, ...}{this function is not intended to be called directly.}
}
\details{
This function produces a 3-d contour or perspective plot. It is intended
to be called by \code{\link{plot.polymars}}.}
\references{Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Martin O'Connor.}
\seealso{\code{\link{polymars}}, 
\code{\link{plot.polymars}}.}
\keyword{smooth}
\keyword{nonlinear}
