\name{lhsc-package}
\alias{lhsc-package}
\docType{package}
\title{
Leaky Hockey Stick Classifier
}
\description{
Leaky hockey stick classifier (LHSC) is a classification algorithm.
This package provides procedures for solving LHSC in reproducing kernel Hilbert spaces. The algorithm is based on the majorization-minimization (MM) principle to compute the entire solution path at a given fine grid of regularization parameters.\cr
}

\details{
Suppose \code{x} is predictor and \code{y} is a binary response. The package computes the entire solution path over a grid of \code{lambda} values.

The main functions of the package \code{\link{lhsc}} include:\cr
\code{lhsc}\cr
\code{cv.lhsc}\cr
\code{predict.lhsc}\cr
\code{plot.lhsc}\cr
\code{plot.cv.lhsc}\cr
}

\author{Oh-Ran Kwon and Hui Zou\cr
Maintainer: Oh-Ran Kwon  \email{kwon0085@umn.edu}}
\references{
Kwon, O. and Zou, H. (2023+)
``Leaky Hockey Stick Loss: The First Negatively Divergent Margin-based Loss Function for Classification"\cr
}

\keyword{lhsc}
\keyword{kernel methods}
\keyword{classification}
\keyword{package}

