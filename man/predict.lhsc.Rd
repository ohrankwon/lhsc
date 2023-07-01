\name{predict.lhsc}
\alias{predict.lhsc}
\title{predict class labels for new observations}
\description{
Predict the binary class labels or the fitted values of an \code{\link{lhsc}} object.
}

\usage{
\method{predict}{lhsc}(object, kern, x, newx, type=c("class", "link"), ...)
}

\arguments{
		\item{object}{A fitted \code{\link{lhsc}} object.}
    \item{kern}{The kernel function used when fitting the \code{\link{lhsc}} object.}
    \item{x}{The predictor matrix, i.e., the \code{x} matrix used when fitting the \code{\link{lhsc}} object.}
		\item{newx}{A matrix of new values for \code{x} at which predictions are to be made. We note that \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
		\item{type}{\code{"class"} or \code{"link"}? \code{"class"} produces the predicted binary class labels and \code{"link"} returns the fitted values. Default is \code{"class"}.}
    \item{\dots}{Not used. Other arguments to \code{predict}.} 
}
 
\details{If \code{"type"} is \code{"class"}, the function returns the predicted class labels. If \code{"type"} is \code{"link"}, the result is \eqn{\beta_0 + x_i'\beta} for the linear case and \eqn{\beta_0 + K_i'\alpha} for the kernel case.}

\value{Returns either the predicted class labels or the fitted values, depending on the choice of \code{type}.}

\author{Oh-Ran Kwon and Hui Zou\cr
Maintainer: Oh-Ran Kwon  \email{kwon0085@umn.edu}}
 
\seealso{\code{\link{lhsc}}}
\examples{
data(BUPA)
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)
lambda = 10^(seq(-3, 3, length.out=10))
kern = rbfdot(sigma=sigest(BUPA$X))
m1 = lhsc(BUPA$X, BUPA$y, kern,
  lambda=lambda, eps=1e-5, maxit=1e5)
predict(m1, kern, BUPA$X, tail(BUPA$X))
}

