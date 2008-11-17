\name{logspline}
\alias{logspline}
\title{ Logspline Density Estimation }
\description{Fits a \code{logspline} density using splines to approximate the log-density
using 
the 1997 knot addition and deletion algorithm (\code{\link{logspline}}). 
The 1992 algorithm is available using the \code{\link{oldlogspline}} function.
}
\usage{logspline(x, lbound, ubound, maxknots = 0, knots, nknots = 0, penalty,
silent = TRUE, mind = -1, error.action = 2) }
\arguments{
\item{x}{ data vector. The data needs to be uncensored. \code{\link{oldlogspline}}
can deal with right- left- and interval-censored data.  }
\item{lbound,ubound}{
lower/upper bound for the support of the density. For example, if there 
is a priori knowledge that the density equals zero to the left of 0, 
and has a discontinuity at 0, 
the user could specify \code{lbound = 0}. However, if the density is  
essentially zero near 0, one does not need to specify \code{lbound}. 
}
\item{maxknots}{
the maximum number of knots. The routine stops adding knots 
when this number of knots is reached. 
The method has an automatic rule 
for selecting maxknots if this parameter is not specified. 
}
\item{knots}{
ordered vector of values (that should cover the complete range of the 
observations), which forces the method to start with these knots. 
Overrules knots.	 
If \code{knots} is not specified, a default knot-placement rule is employed. 
}
\item{nknots}{
forces the method to start with \code{nknots} knots. 
The method has an automatic rule 
for selecting \code{nknots} if this parameter is not specified. 
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses 
the number of knots that minimizes
\code{-2 * loglikelihood + penalty * (number of knots - 1)}. 
The default 
is to use a penalty parameter of \code{penalty = log(samplesize)} as in BIC. The effect of 
this parameter is summarized in \code{\link{summary.logspline}}. 
}
\item{silent}{
should diagnostic output be printed? 
}
\item{mind}{
minimum distance, in order statistics, between knots. 
}
\item{error.action}{how should \code{logspline} deal with non-convergence problems? Very-very rarely
in some extreme situations
\code{logspline} has convergence problems. The only two situations that I am aware of are when
there is effectively a sharp bound, but this bound was not specified, or when the data is severly
rounded. \code{logspline} can deal with this in three ways. If \code{error.action} is 2, the same
data is rerun with the slightly more stable, but less flexible \code{oldlogspline}. The object is translated
in a \code{logspline} object using \code{oldlogspline.to.logspline}, so this is almost
invisible to the user. It is particularly useful when you run simulation studies, as he code can
seemlessly continue. Only the \code{lbound} and \code{ubound} options are passed on to
\code{oldlogspline}, other options revert to the default. If \code{error.action} is 1, a warning is printed,
and \code{logspline} returns nothing (but does not crash). This is useful if you run a
simulation, but do not like to revert to \code{oldlogspline}. If \code{error.action} is 0, the
code crashes using the \code{stop} function.}
}
\value{Object of the class \code{logspline}, that is intended as input for
\code{\link{plot.logspline}} (summary plots), 
\code{\link{summary.logspline}} (fitting summary), 
\code{\link{dlogspline}} (densities), 
\code{\link{plogspline}} (probabilities),
\code{\link{qlogspline}} (quantiles),
\code{\link{rlogspline}} (random numbers from the fitted distribution).

The object has the following members: 

\item{call}{the command that was executed.}
\item{nknots}{the number of knots in the model that was selected.}
\item{coef.pol}{coefficients of the polynomial part of the spline.
The first coefficient is the constant term and
the second is the linear term.}
\item{coef.kts}{coefficients of the knots  part of the spline.
 The \code{k}-th element is the coefficient 
of \eqn{(x-t(k))^3_+} (where \eqn{x^3_+} means the positive part of the third power
of \eqn{x}, 
and \eqn{t(k)} means knot \code{k}).}
\item{knots}{vector of the locations of the knots in the \code{logspline} model.}
\item{maxknots}{the largest number of knots minus one considered during fitting
(i.e. with \code{maxknots = 6} the maximum number of knots is 5).}
\item{penalty}{the penalty that was used.  }
\item{bound}{
first element: 0 - \code{lbound} was \eqn{-\inf}{-infinity,} 1 it was something else; second 
element: \code{lbound}, if specified; third element: 0 - \code{ubound} was \eqn{\inf}{infinity}, 
1 it was something else; fourth element: \code{ubound}, if specified. 
}
\item{samples}{the sample size.}
\item{logl}
{matrix with 3 columns. Column one: number of knots; column two:
model fitted during addition (1) or deletion (2); column 3: log-likelihood.}
\item{range}{range of the input data.}
\item{mind}{minimum distance in order statistics between knots required during fitting
(the actual minimum distance may be much larger).}
}
\references{Charles Kooperberg and Charles J. Stone.  Logspline density estimation
for censored data (1992). \emph{Journal of Computational and Graphical
Statistics}, \bold{1}, 301--328.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{ 
\code{\link{plot.logspline}}, 
\code{\link{summary.logspline}},
\code{\link{dlogspline}}, 
\code{\link{plogspline}},
\code{\link{qlogspline}}, \cr
\code{\link{rlogspline}},
\code{\link{oldlogspline},}
\code{\link{oldlogspline.to.logspline}}.}
\examples{
y <- rnorm(100)
fit <- logspline(y)       
plot(fit)
#
# as (4 == length(-2, -1, 0, 1, 2) -1), this forces these initial knots,
# and does no knot selection
fit <- logspline(y, knots = c(-2, -1, 0, 1, 2), maxknots = 4, penalty = 0)  
#
# the following example give one of the rare examples where logspline
# crashes, and this shows the use of error.action = 2.
#
set.seed(118)
zz <- rnorm(300)
zz[151:300] <- zz[151:300]+5
zz <- round(zz)
fit <- logspline(zz)
#
# you could rerun this with 
# fit <- logspline(zz, error.action=0)
# or
# fit <- logspline(zz, error.action=1)
}
\keyword{distribution}
\keyword{smooth}
