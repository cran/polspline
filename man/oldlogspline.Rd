\name{oldlogspline}
\alias{oldlogspline}
\title{ Logspline Density Estimation - 1992 version }
\description{Fits a \code{logspline} density using splines to approximate the log-density
using 
the 1992 knot deletion algorithm (\code{\link{oldlogspline}}). 
The 1997 algorithm using knot
deletion and addition is available using the \code{\link{logspline}} function.
}
\usage{oldlogspline(uncensored, right, left, interval, lbound = -Inf,
ubound = Inf, nknots, knots, penalty, delete = TRUE) }
\arguments{
\item{uncensored}{
vector of uncensored observations from the distribution whose density is 
to be estimated. If there are no uncensored observations, this argument can 
be omitted. However, either \code{uncensored} or \code{interval} must be specified. 
}
\item{right}{
vector of right censored observations from the distribution 
whose density is to be estimated. If there are no right censored 
observations, this argument can be omitted. 
}
\item{left}{
vector of left censored observations from the distribution 
whose density is to be estimated. If there are no left censored 
observations, this argument can be omitted. 
}
\item{interval}{
two column matrix of lower and upper bounds of observations 
that are interval censored from the distribution whose density is 
to be estimated. If there are no interval censored observations, this 
argument can be omitted. 
}
\item{lbound,ubound}{
lower/upper bound for the support of the density. For example, if there 
is a priori knowledge that the density equals zero to the left of 0, 
and has a discontinuity at 0, 
the user could specify \code{lbound = 0}. However, if the density is  
essentially zero near 0, one does not need to specify \code{lbound}. 
}
\item{nknots}{
forces the method to start with nknots knots (\code{delete = TRUE}) or to fit a 
density with nknots knots (\code{delete = FALSE}). The method has an automatic rule 
for selecting nknots if this parameter is not specified. 
}
\item{knots}{
ordered vector of values (that should cover the complete range of the 
observations), which forces the method to start with these knots (\code{delete = TRUE})
or to fit a density with these knots \code{delete = FALSE}). Overrules \code{nknots}. 
If \code{knots} is not specified, a default knot-placement rule is employed. 
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses 
the number of knots that minimizes \code{-2 * loglikelihood + penalty * (number of knots - 1)}. 
The default is to use  a penalty parameter of \code{penalty = log(samplesize)} as in BIC. The effect of 
this parameter is summarized in \code{\link{summary.oldlogspline}}. 
}
\item{delete}{
should stepwise knot deletion be employed? 
}
}
\value{Object of the class \code{oldlogspline}, that is intended as input for
\code{\link{plot.oldlogspline}}, 
\code{\link{summary.oldlogspline}}, 
\code{\link{doldlogspline}} (densities), 
\code{\link{poldlogspline}} (probabilities),\cr
\code{\link{qoldlogspline}} (quantiles),
\code{\link{roldlogspline}} (random numbers from the fitted distribution).

The object has the following members: 

\item{call}{
the command that was executed. 
}
\item{knots }{
vector of the locations of the knots in the \code{oldlogspline} model. 
old
}
\item{coef}{
coefficients of the spline. The first coefficient is the constant term, 
the second is the linear term and the k-th \eqn{(k>2)} is the coefficient 
of \eqn{(x-t(k-2))^3_+} (where \eqn{x^3_+} means the positive part of the third power
of \eqn{x}, 
and \eqn{t(k-2)} means knot \eqn{k-2}). If a coefficient is zero the corresponding 
knot was deleted from the model. 
}
\item{bound}{
first element: 0 - \code{lbound} was \eqn{-\inf}{-infinity,} 1 it was something else; second 
element: \code{lbound}, if specified; third element: 0 - \code{ubound} was \eqn{\inf}{infinity}, 
1 it was something else; fourth element: \code{ubound}, if specified. 
}
\item{logl}{
the \code{k}-th element is the log-likelihood of the fit with \code{k+2} knots. 
}
\item{penalty}{
the penalty that was used. 
}
\item{sample}{
the sample size that was used. 
}
\item{delete}{
was stepwise knot deletion employed? 
}
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
\code{\link{logspline}},      
\code{\link{oldlogspline}},
\code{\link{plot.oldlogspline}},
\code{\link{summary.oldlogspline}},\cr
\code{\link{doldlogspline}},
\code{\link{poldlogspline}},
\code{\link{qoldlogspline}},
\code{\link{roldlogspline}}.}
\examples{
# A simple example
y <- rnorm(100)
fit <- oldlogspline(y)       
plot(fit)
# An example involving censoring and a lower bound
y <- rlnorm(1000)
censoring <- rexp(1000) * 4
delta <- 1 * (y <= censoring)
y[delta == 0] <- censoring[delta == 0]
fit <- oldlogspline(y[delta == 1], y[delta == 0], lbound = 0)
}
\keyword{distribution}
\keyword{smooth}
