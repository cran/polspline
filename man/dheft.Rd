\name{dheft}
\alias{dheft}
\alias{hheft}
\alias{pheft}
\alias{qheft}
\alias{rheft}
\title{Heft: hazard estimation with flexible tails}
\description{Density (\code{dheft}), cumulative probability (\code{pheft}), hazard rate (\code{hheft}),
quantiles (\code{qheft}), and  random samples (\code{rheft}) from a \code{\link{heft}} object}
\usage{dheft(q, fit) 
hheft(q, fit) 
pheft(q, fit) 
qheft(p, fit) 
rheft(n, fit) }
\arguments{
\item{q}{ vector of quantiles. Missing values (\code{NA}s) are allowed.  }
\item{p}{ vector of probabilities. Missing values (\code{NA}s) are allowed.  }
\item{n}{
sample size. If \code{length(n)} is larger than 1, then
\code{length(n)} random values are returned.
}
\item{fit}{
 \code{heft} object, typically obtained from \code{\link{heft}}.  }
}
\value{Densities (\code{dheft}), hazard rates (\code{hheft}),
probabilities (\code{pheft}), quantiles (\code{qheft}),
or a random sample (\code{rheft}) from a \code{\link{heft}} object.}
\details{
Elements of \code{q} or \code{p} that are missing will cause the
corresponding elements of the result to be missing.
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{ \code{\link{heft}}, \code{\link{plot.heft}}, \code{\link{summary.heft}}.  }
\examples{
fit <- heft(testhare[,1],testhare[,2])
dheft(0:10,fit)
hheft(0:10,fit)
pheft(0:10,fit)
qheft((1:19)/20,fit)
rheft(10,fit)
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
