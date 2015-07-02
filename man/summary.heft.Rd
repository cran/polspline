\name{summary.heft}
\alias{summary.heft}
\alias{print.heft}
\title{Heft: hazard estimation with flexible tails}
\description{This function summarizes both the stepwise selection process of the
model fitting by \code{\link{heft}}, as well as the final model
that was selected using AIC/BIC.}
\usage{\method{summary}{heft}(object, ...) 
\method{print}{heft}(x, ...) }
\arguments{
\item{object,x}{ \code{heft} object, typically the result of \code{\link{heft}}.  }
\item{...}{other arguments are ignored.}
}
\details{These function produce identical printed output. The main body 
is a table with six columns: 

the first column is a possible number 
of knots for the fitted model; 

the second column is 0 if the model was fitted during the 
addition stage and 1 if the model was fitted during the deletion stage; 

the third column is the log-likelihood for the fit; 

the fourth column is \code{-2 * loglikelihood + penalty * (dimension)},
which is the AIC criterion - \code{heft} selected the model with
the minimum value of AIC;

the fifth and sixth columns give the 
endpoints of the interval of values of penalty that would yield the 
model with the indicated number of knots. (\code{NA}s imply that the model is 
not optimal for any choice of penalty.)  

At the bottom of the table the 
number of knots corresponding to the selected model is reported, as are 
the value of penalty that was used and the coefficients of the log-based 
terms in the fitted model and their standard errors. 
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{
\code{\link{heft}},
\code{\link{plot.heft}},
\code{\link{dheft}},
\code{\link{hheft}},
\code{\link{pheft}},
\code{\link{qheft}},
\code{\link{rheft}}.}
\examples{
fit1 <- heft(testhare[,1], testhare[,2])
summary(fit1)
# modify tail behavior
fit2 <- heft(testhare[,1], testhare[,2], leftlog = FALSE, rightlog = FALSE, 
    leftlin = TRUE)   
summary(fit2)
fit3 <- heft(testhare[,1], testhare[,2], penalty = 0)   # select largest model
summary(fit3)
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
