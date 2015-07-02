\name{summary.polyclass}
\alias{summary.polyclass}
\alias{print.polyclass}
\title{Polyclass: polychotomous regression and multiple classification}
\description{This function summarizes both the stepwise selection process of the
model fitting by \code{\link{polyclass}}, as well as the final model
that was selected}
\usage{\method{summary}{polyclass}(object, ...) 
\method{print}{polyclass}(x, ...) }
\arguments{
\item{object,x}{  \code{polyclass} object, typically the result of \code{\link{polyclass}}.  }
\item{...}{other arguments are ignored.}
}
\value{These function summarize a \code{polyclass} fit identically. They also give information 
about fits that could have been obtained with other  
model selection options in \code{\link{polyclass}}. 
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
\code{\link{beta.polyclass}},
\code{\link{cpolyclass}},
\code{\link{ppolyclass}},
\code{\link{rpolyclass}}.}
\examples{
data(iris)
fit.iris <- polyclass(iris[,5], iris[,1:4])
summary(fit.iris)
}
\keyword{smooth}
\keyword{nonlinear}
