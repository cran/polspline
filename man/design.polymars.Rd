\name{design.polymars}
\alias{design.polymars}
\title{Polymars: multivariate adaptive polynomial spline regression}
\description{Produces a design matrux for a model of class \code{polymars}.}
\usage{design.polymars(object, x) }
\arguments{
\item{object}{
object of the class \code{polymars}, typically the result of \code{\link{polymars}}.}
\item{x}{
the predictor values at which the design matrix will be computed.  The 
predictor values can be in a number of formats. It can take the form of a 
vector of length equal to the number of predictors in the original data set 
or it can be shortened to the length of only those predictors that occur in 
the model, in the same order as they appear in the original data set.  
Similarly, \code{x} can take the form of a matrix with the number of columns equal to 
the number of predictors in the original data set, or shortened to the 
number of predictors in the model. 
}
}
\value{The design matrix corresponding to the fitted \code{\link{polymars}} model.}
\references{ Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg}
\seealso{
\code{\link{polymars}},
\code{\link{plot.polymars}},
\code{\link{predict.polymars}},
\code{\link{summary.polymars}}.}
\examples{
data(state)
state.pm <- polymars(state.region, state.x77, knots = 15, classify = TRUE, gcv = 1)
desmat <- design.polymars(state.pm, state.x77)
# compute traditional summary of the fit for the first class
summary(lm(((state.region=="Northeast")*1) ~ desmat -1))
}
\keyword{smooth}
\keyword{nonlinear}
