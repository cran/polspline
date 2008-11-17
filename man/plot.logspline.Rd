\name{plot.logspline}
\alias{plot.logspline}
\title{Logspline Density Estimation }
\description{Plots a \code{logspline} density, distribution function, hazard 
function or survival function 
from 
a logspline density that was fitted using
the 1997 knot addition and deletion algorithm (\code{\link{logspline}}). 
The 1992 algorithm is available using the \code{\link{oldlogspline}} function.
}
\usage{\method{plot}{logspline}(x, n = 100, what = "d", add = FALSE, xlim, xlab = "",
ylab = "", type = "l", ...) 
}
\arguments{
\item{x}{\code{logspline} object, typically the result of \code{\link{logspline}}.}
\item{n}{the number of equally spaced points at which to plot the density.  }
\item{what}{what should be plotted: 
\code{"d"} (density), \code{"p"} (distribution function), \code{"s"} (survival 
function) or \code{"h"} (hazard function).  }
\item{add}{should the plot be added to an existing plot.}
\item{xlim}{
range of data on which to plot. Default is from the 1th to the 99th percentile of
the density, extended by 10\% on each end.}
\item{xlab,ylab}{labels plotted on the axes.  }
\item{type}{type of plot.}
\item{...}{other plotting options, as desired}
}
\details{This function produces a plot of a \code{\link{logspline}} fit at \code{n} equally 
spaced points roughly covering the support of the density. (Use 
\code{xlim = c(from, to)} to change the range of these points.) 
}
\references{Charles Kooperberg and Charles J. Stone.  Logspline density estimation
for censored data (1992). \emph{Journal of Computational and Graphical
Statistics}, \bold{1}, 301--328.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{
\code{\link{logspline}},      
\code{\link{summary.logspline}},
\code{\link{dlogspline}},
\code{\link{plogspline}},
\code{\link{qlogspline}},
\code{\link{rlogspline}},

\code{\link{oldlogspline}}.}
\examples{
y <- rnorm(100)
fit <- logspline(y)       
plot(fit) 
}
\keyword{distribution}
\keyword{smooth}
