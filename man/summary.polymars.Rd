\name{summary.polymars}
\alias{summary.polymars}
\alias{print.polymars}
\title{Polymars: multivariate adaptive polynomial spline regression}
\description{ Gives details of a \code{polymars} object.}
\usage{\method{summary}{polymars}(object, ...) 
\method{print}{polymars}(x, ...) }
\arguments{
\item{object,x}{ object of the class \code{polymars}, typically the result of \code{\link{polymars}}.}
\item{...}{other arguments are ignored.}
}
\details{
These two functions provide identical printed information.
about the fitting steps and the model 
selected. The first data frame contains a row for each step of the fitting 
procedure. In the columns are: a 1 for an addition step or a 0 for a deletion 
step, the size of the model at each step, residual sums of squares (RSS) and  
the generalized cross validation value (GCV), testset residual sums of squares 
or testset misclassification, whatever was used for the model selection. 
The second data frame, model, contains a row for each basis function of the 
model. Each row corresponds to one basis 
function (with two possible components). The pred1 column contains the indices 
of the first predictor of the basis function. Column knot1 is a possible knot  
in this predictor. If this column is NA, the first component is linear. If  
any of the basis functions of the model is categorical then there will be a  
level1 column. Column pred2 is the possible second predictor involved (if  
it is NA the basis function only depends on one  
predictor). Column knot2 contains the possible knot for the predictor pred2,  
and it is NA when this component is linear. This is a similar format  
to the startmodel argument together with an additional first row corresponding to the  
intercept but the startmodel doesn't use a separate column to specify levels of a  
categorical variable . If any predictor in pred2 is categorical then there will be a level2  
column. The column "coefs" (more than one column in the case of multiple response  
regression) contains the coefficients. 
}
\references{ Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Martin O'Connor.}
\seealso{
\code{\link{polymars}},
\code{\link{design.polymars}},
\code{\link{persp.polymars}},
\code{\link{plot.polymars}},
\code{\link{predict.polymars}}.}
\examples{
data(state)
state.pm <- polymars(state.region, state.x77, knots = 15, classify = TRUE)
summary(state.pm)
}
\keyword{smooth}
\keyword{nonlinear}
