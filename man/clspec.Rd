\name{clspec}
\alias{clspec}
\alias{dlspec}
\alias{plspec}
\alias{rlspec}
\title{Lspec: logspline estimation of a spectral distribution}
\description{Autocorrelations, autocovariances
(\code{clspec}), spectral densities and line spectrum (\code{dlspec}),
spectral distributions (\code{plspec}) or a
random time series(\code{rlspec}) from a model fitted with \code{\link{lspec}}.}
\usage{clspec(lag, fit, cov = TRUE, mm) 
dlspec(freq, fit) 
plspec(freq, fit, mm) 
rlspec(n, fit, mean = 0, cosmodel = FALSE, mm)}
\arguments{
\item{lag}{ vector of integer-valued lags for which the  autocorrelations or autocorrelations are to be computed.  }
\item{fit}{ \code{lspec} object, typically the result of \code{\link{lspec}}.}
\item{cov}{ compute autocovariances (\code{TRUE}) or autocorrelations (\code{FALSE}).  }
\item{mm}{ number of points used in integration and the fft. Default is the 
smallest power of two larger than \code{max(fit\$sample, max(lag),1024)}  for
\code{clspec} and \code{plspec}
or the smallest power of two larger than \code{max(fit\$sample, n, max(lag),
1024)} for  (\code{rlspec}).  }
\item{freq}{ vector of frequencies. For \code{plspec} frequencies should be between \eqn{-\pi}{-pi} and \eqn{\pi}{pi}.  }
\item{n}{ length of the random time series to be generated.  }
\item{mean}{ mean level of the time series to be generated.  }
\item{cosmodel}{ indicate that the data should be generated from a model with constant 
harmonic terms rather than a true Gaussian time series.  }
}
\value{Autocovariances or autocorrelations (\code{clspec}); 
values of the spectral distribution at the requested frequencies. (\code{plspec}); 
random time series of length \code{n} (\code{rlspec}); 
or a  list with three components  (\code{dlspec}):
\item{d}{the  spectral density evaluated at the vector of frequencies,}
\item{modfreq}{modified frequencies of the form \eqn{\frac{2\pi j}{T}}{(2pi*j/T)} that are close to the 
frequencies that were requested,}
\item{m}{mass of the line spectrum at the modified frequencies.}
}
\references{ 
Charles Kooperberg, Charles J. Stone, and Young K. Truong (1995). 
Logspline Estimation of a Possibly Mixed Spectral Distribution.
\emph{Journal of Time Series Analysis}, \bold{16}, 359-388.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong. 
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{\code{\link{lspec}}, \code{\link{plot.lspec}}, \code{\link{summary.lspec}}.}
\examples{
data(co2)
co2.detrend <- lm(co2~c(1:length(co2)))$residuals
fit <- lspec(co2.detrend)
clspec(0:12,fit)
plspec((0:314)/100, fit)
dlspec((0:314)/100, fit)
rlspec(length(co2),fit)
}
\keyword{ts}
\keyword{smooth}
