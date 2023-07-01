\name{cv.lhsc}
\alias{cv.lhsc}
\title{cross-validation}
\description{Carry out a cross-validation for \code{\link{lhsc}} to find optimal values of the tuning parameter \code{lambda}.}
\usage{
cv.lhsc(x, y, kern, lambda, nfolds=5, foldid, ...)
}
\arguments{
    \item{x}{A matrix of predictors, i.e., the matrix \code{x} used in \code{\link{lhsc}}.}
    \item{y}{A vector of binary class labels, i.e., the \code{y} used in \code{\link{lhsc}}. \code{y} has to be two levels.}
    \item{kern}{A kernel function.}
    \item{lambda}{A user specified \code{lambda} candidate sequence for cross-validation.}
    \item{nfolds}{The number of folds. Default value is 5. The allowable range is from 3 to the sample size.}
    \item{foldid}{An optional vector with values between 1 and \code{nfold}, representing the fold indices for each observation. If supplied, \code{nfold} can be missing.}
    \item{\dots}{Other arguments being passed to \code{\link{lhsc}}.}
}

\details{This function computes the mean cross-validation error and the standard error by fitting \code{\link{lhsc}} with every fold excluded alternatively.
}
\value{A \code{\link{cv.lhsc}} object including the cross-validation results is return.
    \item{lambda}{The \code{lambda} sequence used in \code{\link{lhsc}}.}
    \item{cvm}{A vector of length \code{length(lambda)}: mean cross-validated error.}
    \item{cvsd}{A vector of length \code{length(lambda)}: estimates of standard error of \code{cvm}.}
    \item{cvupper}{The upper curve: \code{cvm + cvsd}.}
    \item{cvlower}{The lower curve: \code{cvm - cvsd}.}
    \item{lambda.min}{The \code{lambda} incurring the minimum cross validation error \code{cvm}.}
    \item{lambda.1se}{The largest value of \code{lambda} such that error is within one standard error of the minimum.}
    \item{cvm.min}{The cross-validation error corresponding to \code{lambda.min}, i.e., the least error.}
    \item{cvm.1se}{The cross-validation error corresponding to \code{lambda.1se}.}
}


\author{Oh-Ran Kwon and Hui Zou\cr
Maintainer: Oh-ran Kwon  \email{kwon0085@umn.edu}}
\references{
Kwon, O. and Zou, H. (2023+)
``Leaky Hockey Stick Loss: The First Negatively Divergent Margin-based Loss Function for Classification"\cr
}

\seealso{\code{\link{lhsc}} and \code{\link{plot.cv.lhsc}}}
\examples{
set.seed(1)
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(3, -3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m.cv = cv.lhsc(BUPA$X, BUPA$y, kern, lambda=lambda, eps=1e-5, maxit=1e5)
m.cv$lambda.min
}

