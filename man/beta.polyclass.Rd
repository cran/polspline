\name{beta.polyclass}
\alias{beta.polyclass}
\title{Polyclass: polychotomous regression and multiple classification}
\description{Produces a beta-plot for a \code{polyclass} object.  }
\usage{beta.polyclass(fit, which, xsp = 0.4, cex) }
\arguments{
\item{fit}{\code{polyclass} object, typically the result of \code{\link{polyclass}}.  }
\item{which}{which classes should be compared? Default is to compare all classes.  }
\item{xsp}{location of the vertical line to the left of the axis. Useful for making 
high quality, device dependent, graphics.  }
\item{cex}{character size. Default is whatever the present character size is. 
Useful for making high quality, device dependent, graphics. 
}
}
\value{
A beta plot. One line for each basis function. The left part of the plot 
indicates the basis function, the right half the relative location of the 
betas (coefficients) of that basis function, normalized with respect 
to parent basis functions, for all classes. The scaling is supposed 
to suggest a relative importance of the basis functions. This may 
suggest which basis functions are important for separating particular classes. 
}
\references{ Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997). 
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong. 
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{ 
\code{\link{polyclass}},
\code{\link{plot.polyclass}}, 
\code{\link{summary.polyclass}}, 
\code{\link{cpolyclass}}, 
\code{\link{ppolyclass}}, 
\code{\link{rpolyclass}}.}
\note{This is not a generic function, and the complete name, beta.polyclass, has to
be specified.}
\examples{
data(iris)
fit.iris <- polyclass(iris[,5], iris[,1:4])
beta.polyclass(fit.iris)
}
\keyword{smooth}
\keyword{nonlinear}
