\name{summary.hare}
\alias{summary.hare}
\alias{print.hare}
\title{Hare: hazard regression}
\description{This function summarizes both the stepwise selection process of the
model fitting by \code{\link{hare}}, as well as the final model
that was selected using AIC/BIC.}
\usage{summary.hare(object, ...) 
print.hare(x, ...) }
\arguments{
\item{object,x}{ \code{hare} object, typically the result of \code{\link{hare}}.  }
\item{...}{other arguments are ignored.}
}
\details{These function produce identical printed output. The main body consists of 
two tables. 

The first table has six columns: the first column is a 
possible number of dimensions for the fitted model; 

the second column indicates whether this model was fitted during 
the addition or deletion stage; 

the third column is the log-likelihood for the fit; 

the fourth column is \code{-2 * loglikelihood + penalty * (dimension)}, 
which is the AIC criterion - \code{hare} selected the model with 
the minimum value of AIC;  

the last two columns give the 
endpoints of the interval of values of penalty that would yield the 
model with the indicated number of dimensions  
(\code{NA}s imply that the model is not optimal for any choice of penalty). 

At the bottom of the first table the 
dimension of the selected model is reported, as is 
the value of penalty that was used. 

Each row of the second table summarizes the information about  
a basis function in 
the final model. It shows the variables involved, the knot locations, the 
estimated coefficient and its standard error and Wald statistic (estimate/SE). 
}
\note{
Since the basis functions are selected in an adaptive fashion, typically
most Wald statistics are larger than (the magical) 2. These statistics
should be taken with a grain of salt though, as they are inflated because
of the adaptivity of the model selection.}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\seealso{\code{\link{hare}},
\code{\link{plot.hare}},
\code{\link{dhare}},
\code{\link{hhare}},
\code{\link{phare}},
\code{\link{qhare}},
\code{\link{rhare}}.}
\author{ Charles Kooperberg \email{clk@fhcrc.org}.}
\examples{
fit <- hare(testhare[,1], testhare[,2], testhare[,3:8]) 
summary(fit) 
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
