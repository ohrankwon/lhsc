\name{plot.cv.lhsc}
\alias{plot.cv.lhsc}
\title{plot the cross-validation curve}
\description{
Plot cross-validation error curves with the upper and lower standard deviations versus log \code{lambda} values.}
\usage{
\method{plot}{cv.lhsc}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{A fitted \code{\link{cv.lhsc}} object.}
	\item{sign.lambda}{Against \code{log(lambda)} (default) or its negative if \code{sign.lambda=-1}.}
	\item{\dots}{Other graphical parameters being passed to \code{plot}.}
}
\details{This function plots the cross-validation error curves.
}

\author{Oh-Ran Kwon and Hui Zou\cr
Maintainer: Oh-Ran Kwon  \email{kwon0085@umn.edu}}

\seealso{\code{\link{cv.lhsc}}.}
\examples{
set.seed(1)
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m.cv = cv.lhsc(BUPA$X, BUPA$y, kern,
  lambda=lambda, eps=1e-5, maxit=1e5)
m.cv
}
