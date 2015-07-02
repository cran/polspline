\name{plot.heft}
\alias{plot.heft}
\title{Heft: hazard estimation with flexible tails}
\description{Plots a density, distribution function, hazard
function or survival function for
a \code{heft} object.}
\usage{\method{plot}{heft}(x, n = 100, what = "d", add = FALSE, xlim, xlab, ylab,
type, ...) }
\arguments{
\item{x}{\code{heft} object, typically the result of \code{\link{heft}}.  }
\item{n}{the number of equally spaced points at which to plot the function.  }
\item{what}{
what should be plotted: \code{"d"} (density), \code{"p"} (distribution function), \code{"s"} (survival
function) or \code{"h"} (hazard function).  }
\item{add}{should the plot be added to an existing plot?  }
\item{xlim}{plotting limits; default is from the maximum of 0
and 10\% before the 1st percentile
to  the minimmum of
10\% further than the 99th percentile and the largest observation.}
\item{xlab,ylab}{labels for the axes. The default is no labels.}
\item{type}{plotting type. The default is lines.}
\item{...}{all other plotting options are passed on.  }
}
\details{This function produces a plot of a \code{\link{heft}} fit at \code{n} equally
spaced points roughly covering the support of the density. (Use
\code{xlim=c(from,to)} to change the range of these points.)
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{
\code{\link{heft}},
\code{\link{summary.heft}},
\code{\link{dheft}},
\code{\link{hheft}},
\code{\link{pheft}},
\code{\link{qheft}},
\code{\link{rheft}}.}
\examples{
fit1 <- heft(testhare[,1], testhare[,2])
plot(fit1, what = "h")
# modify tail behavior
fit2 <- heft(testhare[,1], testhare[,2], leftlog = FALSE, rightlog = FALSE, 
    leftlin = TRUE)   
plot(fit2, what = "h", add = TRUE,lty = 2)
fit3 <- heft(testhare[,1], testhare[,2], penalty = 0)   # select largest model
plot(fit3, what = "h", add = TRUE,lty = 3)
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
