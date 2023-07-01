\name{plot.lhsc}
\alias{plot.lhsc}
\title{plot coefficients}
\description{
Plot the solution paths for a fitted \code{\link{lhsc}} object.}
\usage{
\method{plot}{lhsc}(x, color=FALSE, ...)
}
\arguments{
  \item{x}{A fitted ``\code{\link{lhsc}}" model.}
	\item{color}{If \code{TRUE}, plots the curves with rainbow colors; otherwise, with gray colors (default).}
	\item{\dots}{Other graphical parameters to \code{plot}.}
}
\details{
Plots the solution paths as a coefficient profile plot.
}
\author{Oh-Ran Kwon and Hui Zou\cr
Maintainer: Oh-Ran Kwon  \email{kwon0085@umn.edu}}
\seealso{\code{lhsc}, \code{predict.lhsc}, \code{coef.lhsc}, \code{plot.lhsc}, and \code{cv.lhsc}.}
 
\examples{
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m1 = lhsc(BUPA$X, BUPA$y, kern, 
  lambda=lambda, eps=1e-5, maxit=1e5)
plot(m1, color=TRUE)
}

