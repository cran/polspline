\name{dhare}
\alias{dhare}
\alias{hhare}
\alias{phare}
\alias{qhare}
\alias{rhare}
\title{Hare: hazard regression}
\description{Density (\code{dhare}), cumulative probability (\code{phare}), hazard rate (\code{hhare}), quantiles
(\code{qhare}), and  random samples (\code{rhare}) from 
a \code{\link{hare}} object.}
\usage{dhare(q, cov, fit) 
hhare(q, cov, fit) 
phare(q, cov, fit) 
qhare(p, cov, fit) 
rhare(n, cov, fit) }
\arguments{
\item{q}{ vector of quantiles. Missing values (\code{NA}s) are allowed.  }
\item{p}{ vector of probabilities. Missing values (\code{NA}s) are allowed.  }
\item{n}{
sample size. If \code{length(n)} is larger than 1, then 
\code{length(n)} random values are returned. 
}
\item{cov}{
covariates. There are several possibilities. If a vector of length 
\code{fit\$ncov} is provided, these covariates are used for all elements of \code{p} or 
\code{q} or for all random numbers. If a matrix of dimension \code{length(p)},
\code{length(q)}, or \code{n} by \code{fit\$ncov} is provided, the rows of \code{cov} are 
matched with the elements of \code{p} or \code{q} or every row of \code{cov} has its own 
random number. If a matrix of dimension \code{m} times \code{fit\$ncov} is provided, 
while \code{length(p) = 1} or \code{length(q) = 1} or \code{n = 1}, the single element of \code{p} or \code{q} is 
used \code{m} times, or \code{m} random numbers with different sets of covariates are 
generated. 
}
\item{fit}{ \code{hare} object, typically obtained from \code{\link{hare}}.  }
}
\value{Densities (\code{dhare}), hazard rates (\code{hhare}), 
probabilities (\code{phare}), quantiles (\code{qhare}), 
or a random sample (\code{rhare}) from a \code{\link{hare}} object.}
}
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
\seealso{\code{\link{hare}},
\code{\link{plot.hare}},
\code{\link{summary.hare}}.}
\examples{
fit <- hare(testhare[,1], testhare[,2], testhare[,3:8])
dhare(0:10, testhare[117,3:8], fit)
hhare(0:10, testhare[1:11,3:8], fit)
phare(10, testhare[1:25,3:8], fit)
qhare((1:19)/20, testhare[117,3:8], fit)
rhare(10, testhare[117,3:8], fit)
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
