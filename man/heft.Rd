\name{heft}
\alias{heft}
\title{Heft: hazard estimation with flexible tails}
\description{Hazard estimation using cubic splines to 
approximate the log-hazard function and special functions to allow
non-polynomial shapes in both tails. 
}
\usage{heft(data, delta, penalty, knots, leftlin, shift, leftlog,
rightlog, maxknots, mindist, silent = TRUE) }
\arguments{
\item{data}{
vector of observations. Observations may or may not be right censored. All
observations should be nonnegative.
}
\item{delta}{
binary vector with the same length as \code{data}. Elements of \code{data}
for which the corresponding element of \code{delta} is 0 are assumed to be
right censored, elements of \code{data}
for which the corresponding element of \code{delta} is 1 are assumed to be
uncensored. If \code{delta} is missing, all observations are assumed to be uncensored.
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses
the number of knots that minimizes \code{-2 * loglikelihood + penalty * (dimension)}.
The default is to use \code{penalty = log(samplesize)} as in BIC. The effect of
this parameter is summarized in \code{\link{summary.heft}}.
}
\item{knots}{
ordered vector of values, which forces the method to start with these knots. 
If \code{knots} is not specified, a default knot-placement rule is employed. 
}
\item{leftlin}{
if \code{leftlin} is \code{TRUE} an extra basis-function, which is linear to the left 
of the first knot, is included in the basis. 
If any of \code{data} is exactly 0, the default of \code{leftlin} is \code{TRUE},
otherwise it is \code{FALSE}.
}
\item{shift}{
parameter for the log terms. Default is \code{quantile(data[delta == 1], .75)}. 
}
\item{leftlog}{
coefficient of \eqn{\log \frac x{x + \mbox{shift}}}{\code{log(x/(x+shift))}}, which must be greater than 
\code{-1}. 
(In particular, if \code{leftlog} equals zero no 
\eqn{\log \frac x{x + \mbox{shift}}}{\code{log(x/(x+shift))}} term 
is included.) If \code{leftlog} is missing its maximum likelihood estimate is used. 
If any of \code{data} is exactly zero, \code{leftlog} is  
set to zero. 
}
\item{rightlog}{
coefficient of \eqn{\log (x + \mbox{shift})}{\code{log(x+shift)}}, which must be greater than 
\code{-1}.
(In particular, if \code{leftlog} equals zero no  
\eqn{\log (x + \mbox{shift})}{\code{log(x+shift)}} term
is included.)
If \code{rightlog} is missing its maximum likelihood estimate is used. 
}
\item{maxknots}{
maximum number of knots allowed in the model (default is 
\eqn{4*n^{0.2})}{\code{4 * n^0.2}}, where \eqn{n} is the length of
\code{data}.
}
\item{mindist}{ minimum distance in order statistics between knots. The
default is 5.  }
\item{silent}{
suppresses the printing of diagnostic output about knots added or deleted, 
Rao-statistics, Wald-statistics and log-likelihoods. 
}
}
\value{An object of class
\code{heft}, which is organized to serve as input for \code{\link{plot.heft}},
\code{\link{summary.heft}}, \code{\link{dheft}} (density), \code{\link{hheft}}
(hazard rate), \code{\link{pheft}} (probabilities), \code{\link{qheft}}
(quantiles), and \code{\link{rheft}} (random numbers).
The object is a list with the following members:
\item{knots }{
vector of the locations of the knots in the \code{heft} model. 
}
\item{logl}{
the \code{k}-th element is the log-likelihood of the fit with \code{k} knots. 
}
\item{thetak}{
coefficients of the knot part of the 
spline. The k-th coefficient is the coefficient 
of  \eqn{(x-t(k))^3_+}. If a coefficient is zero the corresponding 
knot was considered and then deleted from the model. 
}
\item{thetap}{
coefficients of the polynomial part of the spline. 
The first element is the constant term and 
the second element is the linear term. 
}
\item{thetal}{
coefficients of the logarithmic terms. The first element equals 
\code{leftlog} and the second element equals \code{rightlog}. 
}
\item{penalty}{
the penalty that was used. 
}
\item{shift}{ parameter used in the definition of the log terms.  }
\item{sample}{ the sample size.  }
\item{logse}{ the standard errors of \code{thetal}.  }
\item{max}{ the largest element of data.  }
\item{ad}{ vector indicating whether a model of this dimension 
was not fit (2), fit during the addition stage (0) or during 
the deletion stage (1).  }
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\seealso{
\code{\link{hare}},
\code{\link{plot.heft}},
\code{\link{summary.heft}},
\code{\link{dheft}},
\code{\link{hheft}},
\code{\link{pheft}},
\code{\link{qheft}},
\code{\link{rheft}}.}
\examples{
fit1 <- heft(testhare[,1], testhare[,2])
# modify tail behavior
fit2 <- heft(testhare[,1], testhare[,2], leftlog = FALSE, rightlog = FALSE, 
          leftlin = TRUE)   
fit3 <- heft(testhare[,1], testhare[,2], penalty = 0)   # select largest model
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
