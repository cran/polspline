\name{polyclass}
\alias{polyclass}
\title{Polyclass: polychotomous regression and multiple classification}
\description{Fit a polychotomous regression and multiple classification 
using linear splines and selected tensor products.  }
\usage{polyclass(data, cov, weight, penalty, maxdim, exclude, include,
additive = FALSE, linear, delete = 2, fit,  silent = TRUE, 
normweight = TRUE, tdata, tcov, tweight, cv, select, loss, seed) }
\arguments{
\item{data}{vector of classes:
\code{data} should ranges over consecutive integers with 0 or 1 as the minimum value. 
}
\item{cov}{covariates: matrix with as many rows as the length of \code{data}.  }
\item{weight}{optional vector of case-weights.  Should have the same length as 
\code{data}.}
\item{penalty}{
the parameter to be used in the AIC criterion if the 
model selection is carried out by AIC.  The program chooses 
the number of knots that minimizes \code{-2 * loglikelihood + penalty * (dimension)}. 
The default is to use \code{penalty = log(length(data))} as in BIC. If the model 
selection is carried out by cross-validation or using a test set, the 
program uses the number of knots that minimizes 
\code{loss + penalty * dimension * (loss for smallest model)}. In this case 
the default of \code{penalty} is 0. 
}
\item{maxdim}{
maximum dimension (default is 
\eqn{\min(n, 4 * n^{1/3}*(cl-1)}{\code{min(n, 4 * n^(1/3) * (cl - 1)}}, where 
\eqn{n} is \code{length(data)} and
\eqn{cl} the number of classes.
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
\item{additive}{should the model selection be restricted to additive models?  }
\item{linear}{
vector indicating for which of the variables no knots should 
be entered. For example, if \code{linear = c(2, 3)} no knots for either covariate 
2 or 3 are entered. 0 represents time. 
}
\item{delete}{
should complete basis functions be deleted at once (2), should 
only individual dimensions be deleted (1) or should only the addition 
stage of the model selection be carried out (0)? 
}
\item{fit}{
\code{polyclass} object. If  \code{fit} is specified, \code{\link{polyclass}} adds 
basis functions starting with those in \code{fit}. 
}
\item{silent}{
suppresses the printing of diagnostic output about basis functions added 
or deleted, Rao-statistics, Wald-statistics and log-likelihoods. 
}
\item{normweight}{
should the weights be normalized so that they average to one? This option 
has only an effect if the model is selected using AIC. 
}
\item{tdata,tcov,tweight}{
test set. Should satisfy the same requirements as \code{data}, \code{cov} and
\code{weight}. If 
all test set weights are one, \code{tweight} can be omitted. If \code{tdata} and \code{tcov} are 
specified, the model selection is carried out using this test set, irrespective 
of the input for \code{penalty} or \code{cv}. 
}
\item{cv}{
in how many subsets should the data be divided for cross-validation? If \code{cv} is 
specified and tdata is omitted, the model selection is carried out by 
cross-validation.
}
\item{select}{
if a test set is provided, or if the model is selected using cross validation, 
should the model be select that minimizes (misclassification) loss (0), that 
maximizes test set log-likelihood (1) or that minimizes test set 
squared error loss (2)? 
}
\item{loss}{
a rectangular matrix specifying the loss function, whose 
size is the number of 
classes times number of actions. 
Used for cross-validation and test set model 
selection. \code{loss[i, j]} contains the loss for 
assigning action \code{j}  to an object whose true class is \code{i}. 
The default is 1 minus the identity matrix. 
\code{loss} does not need to be square.
}
\item{seed}{
optional 
seed for the random number generator that determines the sequence of the 
cases for cross-validation. If the seed has length 12 or more, 
the first twelve elements are assumed to be \code{.Random.seed}, otherwise 
the function \code{\link{set.seed}} is used.  
If \code{seed} is 0 or \code{rep(0, 12)}, it is assumed that the user 
has already provided a (random) ordering. 
If \code{seed} is not provided, while a fit with 
an element \code{fit\$seed} is provided, 
\code{.Random.seed} is set using \code{set.seed(fit\$seed)}. Otherwise 
the present value of \code{.Random.seed} is used. 
}
}
\value{
The output is an object of class \code{polyclass}, organized
to serve as input for \code{\link{plot.polyclass}}, 
\code{\link{beta.polyclass}},
\code{\link{summary.polyclass}}, \code{\link{ppolyclass}} (fitted probabilities),
\code{\link{cpolyclass}} (fitted classes) and \code{\link{rpolyclass}} (random classes). 
The function returns a list with the following members: 

\item{call}{the command that was executed.  }
\item{ncov}{number of covariates.  }
\item{ndim}{number of dimensions of the fitted model.  }
\item{nclass}{number of classes.  }
\item{nbas}{number of basis functions.  }
\item{naction}{number of possible actions that are considered.  }
\item{fcts}{matrix of size \code{nbas x (nclass + 4)}. each row is a basis function. 
First element: first covariate involved (\code{NA} = constant); 

second element: which knot (\code{NA} means: constant or linear); 

third element: second covariate involved (\code{NA} means: this is a function 
of one variable); 

fourth element: knot involved (if the third element is \code{NA}, of no relevance); 

fifth, sixth,...  element: beta (coefficient) for class one, two, ... 
}
\item{knots}{
a matrix with \code{ncov} rows.
Covariate \code{i} has row \code{i+1}, time has row 1. 
First column: number of knots in this dimension;
other columns: the knots, appended with \code{NA}s to make it a matrix.
}
\item{cv}{
in how many sets was the data divided for cross-validation. 
Only provided if \code{method = 2}. 
}
\item{loss}{
the loss matrix used in cross-validation and test set. 
Only provided if \code{method = 1} or \code{method = 2}.
}
\item{penalty}{
the parameter used in the AIC criterion. Only provided if \code{method = 0}.
}
\item{method}{
0 = AIC, 1 = test set, 2 = cross-validation. 
}
\item{ranges}{
column \code{i} gives the range of the \code{i}-th covariate. 
}
\item{logl}{
matrix with eight or eleven columns. Summarizes fits. 
Column one indicates the dimension, column 
column two the AIC or loss value, whichever was 
used during the model selection 
appropriate, column three four and five give the training set log-likelihood, 
(misclassification) loss and squared error loss, columns six to 
eight give the same information for the test set, column nine (or column 
six if \code{method = 0} or \code{method = 2}) indicates whether the 
model was fitted during the addition stage (1) or during the deletion stage (0), 
column ten and eleven (or seven and eight) the minimum and maximum 
penalty parameter for which AIC would have selected this model. 
}
\item{sample}{sample size. 
}
\item{tsample}{the sample size of the test set. Only prvided if \code{method = 1}.  }
\item{wgtsum}{sum of the case weights.  }
\item{covnames}{names of the covariates.  }
\item{classnames}{(numerical) names of the classes.  }
\item{cv.aic}{the penalty value that was determined optimal by 
by cross validation. Only provided if \code{method = 2}.  }
\item{cv.tab}{table with three columns. Column one and two indicate the penalty parameter 
range for which the cv-loss in column three would be realized.  
Only provided if \code{method = 2}.}
\item{seed}{the random seed that was used to determine the order 
of the cases for cross-validation. 
Only provided if \code{method = 2}.}
\item{delete}{
were complete basis functions deleted at once (2), were 
only individual dimensions deleted (1) or was only the addition 
stage of the model selection carried out (0)? 
}
\item{beta}{
moments of basisfunctions. Needed for \code{\link{beta.polyclass}}.
}
\item{select}{
if a test set is provided, or if the model is selected using cross validation, 
was the model selected that minimized (misclassification) loss (0), that 
maximized test set log-likelihood (1) or that minimized test set 
squared error loss (2)? 
}
\item{anova}{
matrix with three columns. The first two elements in a line 
indicate the subspace to which the line refers. The third element indicates 
the percentage of variance explained by that subspace. 
}
\item{twgtsum}{
sum of the test set case weights (only if \code{method = 1}). 
}
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
\code{\link{polymars}},
\code{\link{plot.polyclass}},
\code{\link{summary.polyclass}},
\code{\link{beta.polyclass}},
\code{\link{cpolyclass}},\cr
\code{\link{ppolyclass}},
\code{\link{rpolyclass}}.}
\examples{
data(iris)
fit.iris <- polyclass(iris[,5], iris[,1:4])
}
\keyword{smooth}
\keyword{nonlinear}
