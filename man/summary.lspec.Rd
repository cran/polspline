\name{summary.lspec}
\alias{summary.lspec}
\alias{print.lspec}
\title{Lspec: logspline estimation of a spectral distribution}
\description{Summary of a model fitted with \code{\link{lspec}} }
\usage{\method{summary}{lspec}(object, ...) 
\method{print}{lspec}(x, ...) }
\arguments{ \item{object,x}{ \code{lspec} object, typically the result of \code{\link{lspec}}.}
\item{...}{other options are ignored.}}
\details{These function produce an identical printed summary of an \code{lspec} object.}
\references{
Charles Kooperberg, Charles J. Stone, and Young K. Truong (1995).
Logspline Estimation of a Possibly Mixed Spectral Distribution.
\emph{Journal of Time Series Analysis}, \bold{16}, 359-388.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{\link{lspec}, \link{plot.lspec}, \link{clspec}, \link{dlspec},
\link{plspec}, \link{rlspec}.}
\examples{
data(co2)
co2.detrend <- lm(co2~c(1:length(co2)))$residuals
fit <- lspec(co2.detrend)
summary(fit)
}
\keyword{ts}
\keyword{smooth}
