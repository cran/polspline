\name{cpolyclass}
\alias{cpolyclass}
\alias{ppolyclass}
\alias{rpolyclass}
\title{Polyclass: polychotomous regression and multiple classification}
\description{ Classify new cases (\code{cpolyclass}), compute class probabilities
for new cases (\code{ppolyclass}), and generate random multinomials for new cases
(\code{rpolyclass}) for a \code{\link{polyclass}} model.}
\usage{cpolyclass(cov, fit)
ppolyclass(data, cov, fit) 
rpolyclass(n, cov, fit) }
\arguments{
\item{cov}{ covariates. Should be a matrix with \code{fit\$ncov} columns.  
For \code{rpolyclass} \code{cov} should either have one row, in 
which case all random numbers are based on the same covariates, or \code{n} 
rows in which case each random number has its own covariates.  }
\item{fit}{ \code{polyclass} object, typically the result of \code{\link{polyclass}}. }
\item{data}{
there are several possibilities. If data is a vector with as many elements 
as cov has rows, each element of data corresponds to a row of cov; if 
only one value is given, the probability of being in that class is computed 
for all sets of covariates. If data is omitted, all class probabilities are 
provided. 
}
\item{n}{ number of pseudo random numbers to be generated.  }
}
\value{ Most likely classes (\code{cpolyclass}),
probabilities (\code{cpolyclass}), or
random classes according to the estimated probabilities (\code{rpolyclass}).
}
\references{ Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{
\code{\link{polyclass}},
\code{\link{plot.polyclass}}, 
\code{\link{summary.polyclass}},
\code{\link{beta.polyclass}}.}   
\examples{
data(iris)
fit.iris <- polyclass(iris[,5], iris[,1:4])
class.iris <- cpolyclass(iris[,1:4], fit.iris)
table(class.iris, iris[,5])
prob.setosa <- ppolyclass(1, iris[,1:4], fit.iris)
prob.correct <- ppolyclass(iris[,5], iris[,1:4], fit.iris) 
rpolyclass(100, iris[64,1:4], fit.iris)
}
\keyword{smooth}
\keyword{nonlinear}
