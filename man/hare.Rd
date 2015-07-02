\name{hare}
\alias{hare}
\title{Hare: hazard regression}
\description{Fit a hazard regression model: linear splines are used to model
the baseline hazard, covariates, and interactions. Fitted models
can be, but do not need to be, proportional hazards models.
}
\usage{hare(data, delta, cov, penalty, maxdim, exclude, include, prophaz = FALSE,
additive = FALSE, linear, fit, silent = TRUE) 
}
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
\item{cov}{
covariates: matrix with as many rows as the length of \code{data}. May be omitted 
if there are no covariates. (If there are no covariates, however,
\code{\link{heft}} will provide a more flexible model using cubic splines.)
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses 
the number of knots that minimizes \code{-2 * loglikelihood + penalty * (dimension)}. 
The default is to use \code{penalty = log(samplesize)} as in BIC. The effect of 
this parameter is summarized in \code{\link{summary.hare}}. 
}
\item{maxdim}{
maximum dimension (default is \eqn{6*\mbox{length(data)}^0.2)}{\code{6 * length(data)^0.2}}. 
}
\item{exclude}{
combinations to be excluded - this should be a matrix with 2 
columns - if for example \code{exclude[1, 1] = 2} and \code{exclude[1, 2] = 3} no 
interaction between covariate 2 and 3 is included. 0 represents time. 
}
\item{include}{
those combinations that can be included. Should have the same format 
as \code{exclude}. Only one of \code{exclude} and \code{include} can be specified .
}
\item{prophaz}{ should the model selection be restricted to proportional hazards models?  }
\item{additive}{ should the model selection be restricted to additive models?  }
\item{linear}{
vector indicating for which of the variables no knots should 
be entered. For example, if \code{linear = c(2, 3)} no knots for either covariate 
2 or 3 are entered. 0 represents time. The default is none.
}
\item{fit}{
\code{\link{hare}} object. If  \code{fit} is specified, \code{\link{hare}} adds 
basis functions starting with those in \code{fit}. 
}
\item{silent}{
suppresses the printing of diagnostic output about basis functions added 
or deleted, Rao-statistics, Wald-statistics and log-likelihoods. 
}
}
\value{An object of class
\code{hare}, which is organized to serve as input for \code{\link{plot.hare}}, 
\code{\link{summary.hare}}, \code{\link{dhare}} (conditional density), \code{\link{hhare}}
(conditional hazard rate), \code{\link{phare}} (conditional probabilities), \code{\link{qhare}}
(conditional quantiles), and \code{\link{rhare}} (random numbers). 
The object is a list with the following members: 
\item{ncov}{ number of covariates.  }
\item{ndim}{ number of dimensions of the fitted model.  }
\item{fcts}{ matrix of size \code{ndim x 6}. each row is a basis function. 
First element: first covariate involved (0 means time); 

second element: which knot (0 means: constant (time) or linear (covariate)); 

third element: second covariate involved (\code{NA} means: this is a function of one variable); 

fourth element: knot involved (if the third element is \code{NA}, of no relevance); 

fifth element: beta; 

sixth element: standard error of beta.}
\item{knots}{
a matrix with \code{ncov} rows.
Covariate \code{i} has row \code{i+1}, time has row 1. 
First column: number of knots in this dimension;
other columns: the knots, appended with \code{NA}s to make it a matrix.
}
\item{penalty}{ the parameter used in the AIC criterion.}
\item{max}{ maximum element of survival data.}
\item{ranges}{ column \code{i} gives the range of the \code{i}-th covariate.}
\item{logl}{
matrix with two columns. The \code{i}-th element of the first column  
is the loglikelihood of the model of dimension \code{i}. The second column indicates whether this 
model was fitted during the addition stage (1) or during the deletion stage (0). 
}
\item{sample}{ sample size.}
}
\references{Charles Kooperberg, Charles J. Stone and Young K. Truong (1995).
Hazard regression.  \emph{Journal of the American Statistical
Association}, \bold{90}, 78-94.

Charles J. Stone, Mark Hansen, Charles Kooperberg, and Young K. Truong.
The use of polynomial splines and their tensor products in extended
linear modeling (with discussion) (1997).  \emph{Annals of Statistics},
\bold{25}, 1371--1470.}
\author{ Charles Kooperberg \email{clk@fredhutch.org}.}
\seealso{\code{\link{heft}},
\code{\link{plot.hare}},
\code{\link{summary.hare}},
\code{\link{dhare}},
\code{\link{hhare}},
\code{\link{phare}},
\code{\link{qhare}},
\code{\link{rhare}}.}
\examples{
fit <- hare(testhare[,1], testhare[,2], testhare[,3:8]) 
}
\keyword{distribution}
\keyword{smooth}
\keyword{survival}
