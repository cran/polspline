\name{plot.lspec}
\alias{plot.lspec}
\title{Lspec: logspline estimation of a spectral distribution}
\description{Plots a spectral density function,
line spectrum, or spectral distribution from a model fitted with \code{\link{lspec}} }
\usage{\method{plot}{lspec}(x, what = "b", n, add = FALSE, xlim, ylim, xlab = "", ylab = "",
type, ...) }
\arguments{
\item{x}{\code{lspec} object, typically the result of \code{\link{lspec}}.}
\item{what}{what should be plotted: b (spectral density and line spectrum 
superimposed), d (spectral density function), 
l (line spectrum)  or p  (spectral distribution function).  }
\item{n}{the number of equally spaced points at which to plot the fit;  default is \code{max(100,fit\$sample)}.  }
\item{add}{indicate that the plot should be added to an existing plot.  }
\item{xlim}{X-axis plotting limits: default is \eqn{c(0,\pi)}{c(0,pi)}, except
when what = "p", when the default is \eqn{c(-\pi,\pi)}{c(-pi,pi)}.}
\item{ylim}{Y-axis plotting limits.}
\item{xlab,ylab}{axis labels.}
\item{type}{plotting type; default is \code{"l"} when \code{what = "d"}
 and \code{what = "p"}, \code{"h"} when \code{what = "l"}, and
a combination of \code{"h"} and \code{"l"} when \code{what ="b"}}
\item{...}{all regular plotting options are passed on.  }
}
\note{If \code{what = "p"} 
the plotting range cannot extend beyond the interval \eqn{[-\pi,\pi]}{[-pi,pi]}.  }
\references{
Charles Kooperberg, Charles J. Stone, and Young K. Truong (1995).
Logspline Estimation of a Possibly Mixed Spectral Distribution.
\emph{Journal of Time Series Analysis}, \bold{16}, 359-388.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{\code{\link{lspec}}, \code{\link{summary.lspec}}, \code{\link{clspec}}, \code{\link{dlspec}},
\code{\link{plspec}}, \code{\link{rlspec}}.}
\examples{
data(co2)
co2.detrend <- lm(co2~c(1:length(co2)))$residuals
fit <- lspec(co2.detrend)
plot(fit)
}
\keyword{ts}
\keyword{smooth}
