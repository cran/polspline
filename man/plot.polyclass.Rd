\name{plot.polyclass}
\alias{plot.polyclass}
\title{Polyclass: polychotomous regression and multiple classification}
\description{Probability or classification plots for a  \code{\link{polyclass}} model.  }
\usage{\method{plot}{polyclass}(x, cov, which, lims, what, data, n, xlab="", ylab="",
zlab="", ...) }
\arguments{
\item{x}{\code{polyclass} object, typically the result of \code{\link{polyclass}}. }
\item{cov}{
a vector of length \code{fit\$ncov}, indicating for which combination of  
covariates the plot should be made. Can never be omitted. Should always have 
length \code{fit\$ncov}, even if some values are irrelevant.  }
\item{which}{
for which covariates should the plot be made.  
Number or a character string defining the name, if the 
same names were used with the call to \code{\link{polyclass}}. Which should have length one if 
\code{what} is 6 or larger and length two if \code{what} is 5 or smaller.  }
\item{lims}{plotting limits. If omitted, the plot is made over the same range 
of the covariate as in the original data. Otherwise a vector of 
length two of the form \code{c(min, max)} if what is 6 or larger and a vector of 
length four of the form \code{c(xmin, xmax, ymin ,ymax)} if \code{what} is 5 or smaller.  }
\item{what}{
an integer between 1 and 8, defining the type of plot to be made. 
\enumerate{
\item Plots the probability of one class as a contour plot of two variables. 
\item Plots the probability of one class as a perspective plot of two variables. 
\item Plots the probability of one class as an image plot of two variables. 
\item Classifies the area as a contour plot of two variables. 
\item Classifies the area as an image plot of two variables. 
\item Classifies the line as a plot of one variable. 
\item Plots the probabilities of all classes as a function of one variable. 
\item Plots the probability of one class as a function of one variable. 
}
}
\item{data}{Class for which the plot is made. Should be provided if \code{what} is 1, 2, 3 or 8.  }
\item{n}{
the number of equally spaced points at which to plot the fit. The 
default is 250 if \code{what} is 6 or larger or 50 (which results in 2500 plotting 
points) if \code{what} is 5 or smaller. 
}
\item{xlab,ylab,zlab}{axis plotting labels.}
\item{...}{
all other options are passed on.  }
}
\references{Charles Kooperberg, Smarajit Bose, and  Charles J. Stone (1997).
Polychotomous regression. \emph{Journal of the American Statistical
Association}, \bold{92}, 117--127.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{
\code{\link{polyclass}},
\code{\link{summary.polyclass}},
\code{\link{beta.polyclass}},
\code{\link{cpolyclass}},
\code{\link{ppolyclass}},
\code{\link{rpolyclass}}.}
\examples{
data(iris)
fit.iris <- polyclass(iris[,5], iris[,1:4])
plot(fit.iris, iris[64,1:4], which=c(3,4), data=2, what=1) 
plot(fit.iris,iris[64,1:4], which=c(3,4), what=5) 
plot(fit.iris,iris[64,1:4], which=4, what=7) 
}
\keyword{smooth}
\keyword{nonlinear}
