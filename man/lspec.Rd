\name{lspec}
\alias{lspec}
\title{Lspec: logspline estimation of a spectral distribution}
\description{ Fit an \code{lspec} model
to a time-series or a periodogram.}
\usage{lspec(data, period, penalty, minmass, knots, maxknots, atoms, maxatoms,
maxdim , odd = FALSE, updown = 3, silent = TRUE) }
\arguments{
\item{data}{
time series (exactly one of \code{data} and \code{period} should be specified). 
If \code{data} is specified, \code{lspec} first computes the modulus 
of the fast Fourier transform 
of the series using the function
\code{\link{fft}}, resulting in a periodogram of length
\code{floor(length(data)/2)}. 
}
\item{period}{
value of the periodogram for a time series at frequencies 
\eqn{\frac{2\pi j}T}{2pi*j/T}, for \eqn{1\leq j \leq T/2}{1<=j<=T/2}. If period is specified, odd should indicate 
whether the length of the series T is odd \code{(odd = TRUE)} or even \code{(odd = FALSE)}. 
Exactly one of \code{data} and \code{period} should be specified. 
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses 
the number of basis 
functions that minimizes \code{-2 * loglikelihood + penalty * (number of basis 
functions)}. 
Default is to use a penalty parameter of \code{penalty = log(length(period))} as in BIC.  }
\item{minmass}{
threshold value for atoms. No atoms having smaller mass than \code{minmass} are 
included in the model. If \code{minmass} takes its default value, in 
95\% of the samples, when data is Gaussian white noise, the model will not 
contain atoms.  }
\item{knots}{
ordered vector of values, which forces the method to start with these knots. 
If \code{knots} is not specified, the program starts with one knot at zero and then 
employs stepwise addition of knots and atoms. 
}
\item{maxknots}{ maximum number of knots allowed in the model. Does not need to be 
specified, since the program has a default for \code{maxdim} and 
the number of dimensions equals the number of knots plus the number of 
atoms. If \code{maxknots = 1} the fitted spectral density function is 
constant. 
}
\item{atoms}{
ordered vector of values, which forces the method to start with discrete 
components at these frequencies. The values of atoms are rounded 
to the nearest multiple of \eqn{\frac{2\pi}T}{2*pi/T}. 
If atoms is not specified, the program starts with no atoms and then performs 
stepwise addition of knots and atoms. 
}
\item{maxatoms}{
maximum number of discrete components allowed in the model. Does not need to be 
specified, since the program has a default for \code{maxdim} and 
the number of dimensions equals the number of knots plus the number of 
atoms. If \code{maxatoms = 0} a continuous 
spectral distribution is fit. 
}
\item{maxdim}{
maximum number of basis 
functions allowed in the model (default is 
\eqn{\max(15,4\times\mbox{length(period)}^{0.2})}{max(15,4*length(period)^0.2)}). 
}
\item{odd}{
see \code{period}. If \code{period} is not specified, \code{odd} is not relevant. 
}
\item{updown}{
the maximal number of times that \code{lspec} should go through a cycle of stepwise 
addition and stepwise deletion until a stable solution is reached. 
}
\item{silent}{
should printing of information be suppressed?}
}
\value{Object of class \code{lspec}.
The output is organized to serve as input for \code{\link{plot.lspec}}
(summary plots),
\code{\link{summary.lspec}} (summarizes fitting), \code{\link{clspec}} (for
autocorrelations and autocovariances), \code{\link{dlspec}} (for spectral density and line-spectrum,) 
\code{\link{plspec}} (for the spectral distribution), and \code{\link{rlspec}}
(for random time series with the same spectrum).

\item{call}{
the command that was executed. 
}
\item{thetap}{
coefficients of the polynomial part of the spline. 
}
\item{nknots}{
the number of knots that were retained. 
}
\item{knots}{
vector of the locations of the knots in the logspline model. 
Only the knots that were retained are in this vector.  
}
\item{thetak}{
coefficients of the knot part of the 
spline. The k-th coefficient is the coefficient 
of  \eqn{(x-t(k))^3_+}. 
}
\item{natoms}{
the number of atoms that were retained. 
}
\item{atoms}{
vector of the locations of the atoms in the model. 
Only the atoms that were retained are in this vector.  
}
\item{mass}{
The k-th coefficient is the mass at \code{atom[k]}. 
}
\item{logl}{
the log-likelihood of the model. 
}
\item{penalty}{
the penalty that was used. 
}
\item{minmass}{
the minimum mass for an atom that was allowed. 
}
\item{sample}{
the sample size that was used, either computed as \code{length(data)} or 
as \code{(2 * length(period))}  when \code{odd = FALSE} or as
\code{(2 * length(period) + 1)}  when \code{odd = TRUE}. 
}
\item{updown}{
the actual number of times that \code{lspec} went through a cycle of  
stepwise addition and stepwise  deletion  
until a stable solution was reached, or 
minus the number of times that lspec went through a cycle of  
stepwise addition and stepwise  deletion until it decided to quit. 
}
}
\references{
Charles Kooperberg, Charles J. Stone, and Young K. Truong (1995).
Logspline Estimation of a Possibly Mixed Spectral Distribution.
\emph{Journal of Time Series Analysis}, \bold{16}, 359-388.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{\code{\link{plot.lspec}}, \code{\link{summary.lspec}}, \code{\link{clspec}}, \code{\link{dlspec}},
\code{\link{plspec}}, \code{\link{rlspec}}.}
\examples{
data(co2)
co2.detrend <- unstrip(lm(co2~c(1:length(co2)))$residuals)
fit <- lspec(co2.detrend)
}
\keyword{ts}
\keyword{smooth}
