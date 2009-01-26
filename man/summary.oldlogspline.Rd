\name{summary.oldlogspline}
\alias{summary.oldlogspline}
\alias{print.oldlogspline}
\title{ Logspline Density Estimation - 1992 version }
\description{ This function summarizes both the stepwise selection process of the
model fitting by \code{\link{oldlogspline}}, as well as the final model
that was selected using AIC/BIC. A
 \code{logspline} object was fit using
the 1992 knot deletion algorithm (\code{\link{oldlogspline}}). 
The 1997 algorithm using knot
deletion and addition is available using the \code{\link{logspline}} function.
}
\usage{\method{summary}{oldlogspline}(object, ...) 
\method{print}{oldlogspline}(x, ...)}
\arguments{
\item{object,x}{
\code{oldlogspline} object, typically the result of \code{\link{oldlogspline}}
}
\item{...}{ other arguments are ignored.}
}
\details{
These function produces the same printed output. The main body 
is a table with five columns: the first column is a possible number 
of knots for the fitted model; 

the second column is the log-likelihood for the fit; 

the third column is \code{-2 * loglikelihood + penalty * (number of knots - 1)},
which is the AIC criterion; \code{\link{logspline}} selected the model with
the smallest value of AIC;


the fourth and fifth columns give the 
endpoints of the interval of values of penalty that would yield the 
model with the indicated number of knots. (\code{NA}s imply that the model is 
not optimal for any choice of \code{penalty}.) At the bottom of the table the 
number of knots corresponding to the selected model is reported, as is 
the value of penalty that was used. 
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
\code{\link{doldlogspline}},
\code{\link{poldlogspline}},\cr
\code{\link{qoldlogspline}},
\code{\link{roldlogspline}}.}
\examples{
y <- rnorm(100)
fit <- oldlogspline(y)       
summary(fit) 
}
\keyword{distribution}
\keyword{smooth}
