\name{unstrip}
\alias{unstrip}
\title{Reformat data as vector or matrix}
\description{This function tries to convert a date.frame or a matrix to a no-frills matrix without labels,
and a vector or time-series to a no-frills vector without labels.}
\usage{unstrip(x) }
\arguments{ \item{x}{ one- or two-dimensional object.} }
\value{If \code{x} is two-dimensional a  matrix
without names, if \code{x} is one-dimensional  a numerical vector}
\details{Many of the functions for \code{\link{logspline}}, \code{\link{oldlogspline}}, 
\code{\link{lspec}}, \code{\link{polyclass}},
\code{\link{hare}}, \code{\link{heft}}, and \code{\link{polymars}} were written in the ``before data.frame'' era;
\code{unstrip} attempts to keep all these functions useful with more advanced input objects.
In particular, many of these functions call \code{unstrip} before doing anything else.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\examples{
data(co2)
unstrip(co2)
data(iris)
unstrip(iris)
}
\keyword{utilities}
\keyword{classes}
