\name{doldlogspline}
\alias{doldlogspline}
\alias{poldlogspline}
\alias{qoldlogspline}
\alias{roldlogspline}
\title{Logspline Density Estimation - 1992 version }
\description{Probability density function (\code{doldlogspline}), distribution
function (\code{poldlogspline}), quantiles 
(\code{qoldlogspline}), and  random samples (\code{roldlogspline}) from
a logspline density that was fitted using
the 1992 knot deletion algorithm (\code{\link{oldlogspline}}). 
The 1997 algorithm using knot
deletion and addition is available using the \code{\link{logspline}} function.
}
\usage{doldlogspline(q, fit) 
poldlogspline(q, fit) 
qoldlogspline(p, fit) 
roldlogspline(n, fit) }
\arguments{
\item{q}{
vector of quantiles. Missing values (NAs) are allowed. 
}
\item{p}{
vector of probabilities. Missing values (NAs) are allowed. 
}
\item{n}{
sample size. If \code{length(n)} is larger than 1, then 
\code{length(n)} random values are returned.
}
\item{fit}{
\code{oldlogspline} object, typically the result of \code{\link{oldlogspline}}.
}
}
\value{ Densities (\code{doldlogspline}), probabilities (\code{poldlogspline}), quantiles (\code{qoldlogspline}), 
or a random sample (\code{roldlogspline})
from an \code{oldlogspline} density that was fitted using
knot deletion.
}
\details{
Elements of  \code{q} or \code{p} that are missing will cause the 
corresponding elements of the result to be missing. 
}
\references{Charles Kooperberg and Charles J. Stone.  Logspline density estimation
for censored data (1992). \emph{Journal of Computational and Graphical
Statistics}, \bold{1}, 301--328.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{
\code{\link{logspline}},
\code{\link{oldlogspline}},
\code{\link{plot.oldlogspline}},
\code{\link{summary.oldlogspline}}}
\examples{
x <- rnorm(100)
fit <- oldlogspline(x)
qq <- qoldlogspline((1:99)/100, fit)
plot(qnorm((1:99)/100), qq)                  # qq plot of the fitted density
pp <- poldlogspline((-250:250)/100, fit)
plot((-250:250)/100, pp, type = "l")
lines((-250:250)/100, pnorm((-250:250)/100)) # asses the fit of the distribution
dd <- doldlogspline((-250:250)/100, fit)
plot((-250:250)/100, dd, type = "l")
lines((-250:250)/100, dnorm((-250:250)/100)) # asses the fit of the density
rr <- roldlogspline(100, fit)                # random sample from fit
}
\keyword{distribution}
\keyword{smooth}
