\name{oldlogspline.to.logspline}
\alias{oldlogspline.to.logspline}
\title{Logspline Density Estimation - 1992 to 1997 version }
\description{Translates an \code{oldlogspline} object in an
\code{logspline} object. This routine is mostly used in \code{logspline},
as it allows the routine to use \code{oldlogspline} for some situations
where \code{logspline} crashes. The other use is when you have censored data,
and thus have to use \code{oldlogspline} to fit, but wish to use the
auxiliary routines from \code{logspline}.}
\usage{oldlogspline.to.logspline(obj, data) }
\arguments{
\item{obj}{
object of class \code{logspline}
}
\item{data}{
the original data. Used to compute the \code{range} component of the
new object. if \code{data} is not available, the 1/(n+1) and n/(n+1)
quantiles of the fitted distribution are used for \code{range}.
}
}
\value{ object of the class \code{logspline}. The \code{call} component
of the new object is not useful. The \code{delete} component of the old
object is ignored.}
\references{Charles Kooperberg and Charles J. Stone.  Logspline density estimation
for censored data (1992). \emph{Journal of Computational and Graphical
Statistics}, \bold{1}, 301--328.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{
\code{\link{logspline}},
\code{\link{oldlogspline}}.}
\examples{
x <- rnorm(100)
fit.old <- oldlogspline(x)
fit.translate <- oldlogspline.to.logspline(fit.old,x)
fit.new <- logspline(x)
plot(fit.new)
plot(fit.old,add=TRUE,col=2)
#
# should look almost the same, the differences are the
# different fitting routines
#
}
\keyword{distribution}
\keyword{smooth}
