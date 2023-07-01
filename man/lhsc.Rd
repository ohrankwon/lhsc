\name{lhsc}
\alias{lhsc}
\title{solve linear LHSC and Kernel LHSC}
\description{Fit the LHSC in input space and reproducing kernel Hilbert space. The solution path is computed at a grid of values of tuning parameter \code{lambda}.}
\usage{
lhsc(x, y, kern, lambda, eps=1e-05, maxit=1e+05)
}
\arguments{
	\item{x}{A numerical matrix with \eqn{N} rows and \eqn{p} columns for predictors.}
	\item{y}{A vector of length \eqn{N} for binary responses. The element of \code{y} is either -1 or 1.}
	\item{kern}{A kernel function; see \code{\link{dots}}.}
  \item{lambda}{A user supplied \code{lambda} sequence.}
	\item{eps}{The algorithm stops when \eqn{| \beta^{old} - \beta^{new} |} is less than \code{eps}. Default value is \code{1e-5}.}
	\item{maxit}{The maximum of iterations allowed. Default is 1e5.}
}

\details{The leaky hockey stick loss is \eqn{V(u)=1-u}{V(u) = 1 - u} if \eqn{u \le 1}{u <= 1} and \eqn{-\log u}{-log u} if \eqn{u > 1}{u > 1}. The value of \eqn{\lambda}, i.e., \code{lambda}, is user-specified. 

In the linear case (\code{kern} is the inner product and N > p), the \code{\link{lhsc}} fits a linear LHSC by minimizing the L2 penalized leaky hockey stick loss function,
\deqn{L(\beta_0,\beta) := \frac{1}{N}\sum_{i=1}^N V(y_i(\beta_0 + X_i'\beta)) + \lambda \beta' \beta.}{(1/N) * sum_i [V(y_i(\beta_0 + X_i'\beta))] + \lambda \beta' \beta.} 

If a linear LHSC is fitted when N < p, a kernel LHSC with the linear kernel is actually solved. In such case, the coefficient \eqn{\beta}{\beta} can be obtained from \eqn{\beta = X'\alpha.}{\beta = X'\alpha.} 

In the kernel case, the \code{\link{lhsc}} fits a kernel LHSC by minimizing
\deqn{L(\alpha_0,\alpha) := \frac{1}{n}\sum_{i=1}^n V(y_i(\alpha_0 + K_i' \alpha)) + \lambda \alpha' K \alpha,}{(1/n) * sum_i [V(y_i(\alpha_0 + K_i' \alpha))] + \lambda \alpha' K \alpha,}
where \eqn{K}{K} is the kernel matrix and \eqn{K_i}{K_i} is the ith row. 
}

\value{
An object with S3 class \code{\link{lhsc}}.
  \item{alpha}{A matrix of LHSC coefficients at each \code{lambda} value. The dimension is \code{(p+1)*length(lambda)} in the linear case and \code{(N+1)*length(lambda)} in the kernel case.}
  \item{lambda}{The \code{lambda} sequence.}
  \item{npass}{The total number of FISTA iterations for all lambda values. }
  \item{jerr}{Warnings and errors; 0 if none.}
  \item{info}{A list including parameters of the loss function, \code{eps}, \code{maxit}, \code{kern}, and \code{wt} if a weight vector was used.}
  \item{call}{The call that produced this object.}
}

\author{Oh-ran Kwon and Hui Zou\cr
Maintainer: Oh-ran Kwon  \email{kwon0085@umn.edu}}
\references{
Kwon, O. and Zou, H. (2023+)
``Leaky Hockey Stick Loss: The First Negatively Divergent Margin-based Loss Function for Classification"\cr
}

\seealso{\code{\link{predict.lhsc}}, \code{\link{plot.lhsc}}, and \code{\link{cv.lhsc}}.}
\examples{
data(BUPA)
# standardize the predictors
BUPA$X = scale(BUPA$X, center=TRUE, scale=TRUE)

# a grid of tuning parameters
lambda = 10^(seq(3, -3, length.out=10))

# fit a linear LHSC
kern = vanilladot()
DWD_linear = lhsc(BUPA$X, BUPA$y, kern,
  lambda=lambda, eps=1e-5, maxit=1e5)

# fit a kernel LHSC using Gaussian kernel
kern = rbfdot(sigma=1)
DWD_Gaussian = lhsc(BUPA$X, BUPA$y, kern,
  lambda=lambda, eps=1e-5, maxit=1e5)
}
\keyword{lhsc}
\keyword{kernel methods}
\keyword{classification}
