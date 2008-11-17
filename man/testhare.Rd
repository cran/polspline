\name{testhare}
\alias{testhare}
\title{Fake survival data for Hare and Heft}
\description{Fake survival analysis data set for testing \code{\link{hare}} and \code{\link{heft}}}
\usage{testhare}
\format{
A matrix with 2000 lines (observations) and 8 columns. Column 1 is
intended to be the survival time, column 2 the censoring indicator, and
columns 3 through 8 are predictors (covariates).}
\source{I started out with a real data set; then I sampled, transformed and added noise.
Virtually no number is unchanged.}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{ \code{\link{hare}}, \code{\link{heft}}.}
\examples{
harefit <- hare(testhare[,1], testhare[,2], testhare[,3:8]) 
heftfit <- heft(testhare[,1], testhare[,2])
}
\keyword{survival}
\keyword{datasets}
