\name{predict.polymars}
\alias{predict.polymars}
\title{Polymars: multivariate adaptive polynomial spline regression}
\description{Produces fitted values for a model of class \code{polymars}.}
\usage{\method{predict}{polymars}(object, x, classify = FALSE, intercept, ...) }
\arguments{
\item{object}{
object of the class \code{polymars}, typically the result of \code{\link{polymars}}.}
\item{x}{
the predictor values at which the fitted values will be computed.  The 
predictor values can be in a number of formats. It can take the form of a 
vector of length equal to the number of predictors in the original data set 
or it can be shortened to the length of only those predictors that occur in 
the model, in the same order as they appear in the original data set.  
Similarly, \code{x} can take the form of a matrix with the number of columns equal to 
the number of predictors in the original data set, or shortened to the 
number of predictors in the model. 
}
\item{classify}{
if the original call to polymars was for a classification problem and you would 
like the classifications (class predictions), set this option equal to \code{TRUE}. Otherwise the 
function returns a response column for each class (the highest values in each 
row is its class for the case when \code{classify = TRUE}).  }
\item{intercept}{
 Setting intercept equal to \code{FALSE} evaluates the object without intercept. The  
intercept may also be given any numerical value which overrides the fitted  
coefficient from the object. The defualt is  \code{TRUE}.
}
\item{...}{other arguments are ignored.}
}
\value{A matrix of fitted values.
The number of columns in the 
returned matrix equals the number of responses in the original call to \code{\link{polymars}}.}
\references{ Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Martin O'Connor.}
\seealso{
\code{\link{polymars}},
\code{\link{design.polymars}},
\code{\link{plot.polymars}},
\code{\link{summary.polymars}}.}
\examples{
data(state)
state.pm <- polymars(state.region, state.x77, knots = 15, classify = TRUE, gcv = 1)
table(predict(state.pm, x = state.x77, classify = TRUE), state.region)
}
\keyword{smooth}
\keyword{nonlinear}
