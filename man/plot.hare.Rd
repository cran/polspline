\name{plot.hare}
\alias{plot.hare}
\title{Hare: hazard regression}
\description{Plots a density, distribution function, hazard 
function or survival function for
a \code{hare} object.}
\usage{\method{plot}{hare}(x, cov, n = 100, which = 0, what = "d", time, add = FALSE, xlim,
xlab, ylab, type, ...) }
\arguments{
\item{x}{\code{hare} object, typically the result of \code{\link{hare}}.  }
\item{cov}{a vector of length \code{fit\$ncov}, indicating for which combination of  
covariates the plot should be made. Can be omitted only if \code{fit\$ncov} is 0.  }
\item{n}{the number of equally spaced points at which to plot the function.  }
\item{which}{for which coordinate should the plot be made. 0: time; positive value 
i: covariate i. Note that if which is the positive value i, then the 
element corresponding to this covariate must be given in \code{cov} even 
though its actual value is irrelevant.  }
\item{what}{
what should be plotted: \code{"d"} (density), \code{"p"} (distribution function), \code{"s"} (survival 
function) or \code{"h"} (hazard function).  }
\item{time}{if which is not equal to 0, the value of time for which the plot should be made.  }
\item{add}{should the plot be added to an existing plot?  }
\item{xlim}{plotting limits; default is from the maximum of 0 
and 10\% before the 1st percentile
to  the minimmum of
10\% further than the 99th percentile and the largest observation.}
\item{xlab,ylab}{labels for the axes. Per default no labels are printed.}
\item{type}{plotting type. The default is lines.}
\item{...}{all other plotting options are passed on.  }
}
\details{This function produces a plot of a \code{\link{hare}} fit at \code{n} equally 
spaced points roughly covering the support of the density. (Use 
\code{xlim=c(from,to)} to change the range of these points.) 
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{\code{\link{hare}},
\code{\link{summary.hare}},
\code{\link{dhare}},
\code{\link{hhare}},
\code{\link{phare}},
\code{\link{qhare}},
\code{\link{rhare}}.}
\examples{
fit <- hare(testhare[,1], testhare[,2], testhare[,3:8])       
# hazard curve for covariates like case 1 
plot(fit, testhare[1,3:8], what = "h") 
# survival function as a function of covariate 2, for covariates as case 1 at t=3 
plot(fit, testhare[1,3:8], which = 2, what = "s",  time = 3)  
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
