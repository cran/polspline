\name{plot.polymars}
\alias{plot.polymars}
\title{Polymars: multivariate adaptive polynomial spline regression}
\description{Produces two and three dimensional plots of the 
fitted values from a \code{polymars} object. 
}
\usage{plot.polymars(x, predictor1, response, predictor2, xx, add = FALSE, n,
xyz = FALSE, contour.polymars = FALSE, xlim, ylim, intercept, ...) }
\arguments{
\item{x}{\code{polymars} object, typically the result of \code{\link{polymars}}. }
\item{predictor1}{the index of a predictor that was used when the \code{polymars} model was fit. 
For the two dimensional plots, this variable is plotted along the X-axis.  }
\item{response}{
if the model was fitted to multiple response data the response index 
should be specified. }
\item{predictor2}{
the index of a predictor that was used when the \code{polymars} model was fit. For the three  
dimensional plots, this variable is plotted along the Y-axis.   
See \code{xyz}. 
}
\item{xx}{
should be a vector of length equal to the number of predictors in the  
original data set. The values should be in the same order as in the original  
dataset. By default the function uses the median values of the data that was  
used to fit the model.  Although the values for predictor and predictor2 are  
not used, they should still be provided as part of \code{xx}. 
}
\item{add}{
should the plot be added to a previously created plot? Works only for two 
dimensional plots. 
}
\item{n}{
number of plotting points (2 dimensional plot) or plotting points along each 
axis (3 dimensional plot). The default is \code{n = 100} for 2 dimensional plots and 
\code{n = 33} for 3 dimensional plots. 
}
\item{xyz}{
is the plot being made a 3 dimensional plot? 
If there is only one response it need not be set, if two numerical values  
accompany the model in the call they will be understood as two predictors 
for a 3-d plot. By default a 3-d plot uses the \code{\link{persp}} function. 
Categorical predictors cannot be used for 3 dimensional plots. 
}
\item{contour.polymars}{
if the plot being made a 3 dimensional plot should it be made as a contour plot 
(\code{TRUE}) or a perspective plot (\code{FALSE}).
function \link{contour} is being made. 
}
\item{intercept}{
 Setting intercept equal to \code{FALSE} evaluates the object without intercept. The  
intercept may also be given any numerical value which overrides the fitted  
coefficient from the object. The default is \code{TRUE}.
}
\item{xlim,ylim}{Plotting limits. The function tries to choose intelligent limits itself}
\item{...}{
other options are passed on.
}
}
\details{
This function produces a 2-d plot of 1 predictor and response of a \code{polymars} object
at n equally spaced points or a 3-d plot of two predictors and response of a 
\code{polymars} object.  The range of the plot is by default equal to the range of the 
particular predictor(s) in the original data, but this can be changed by  
\code{xlim = c(from, to)} and
\code{ylim = c(from, to)}.
}
\references{Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Martin O'Connor.}
\seealso{
\code{\link{design.polymars}},
\code{\link{polymars}},
\code{\link{predict.polymars}},
\code{\link{summary.polymars}}.}
\examples{
data(state)
state.pm <- polymars(state.region, state.x77, knots = 15, classify = TRUE, gcv = 1)
plot(state.pm, 3, 4)
}
\keyword{smooth}
\keyword{nonlinear}
