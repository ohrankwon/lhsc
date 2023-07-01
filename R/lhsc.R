lhsc = function(x, y, kern, lambda, eps=1e-05, maxit=1e+05) {
  ####################################################################
  ## data setup
  this.call = match.call()
  if (length(levels(factor(y))) == 2)
    y = c(-1, 1)[as.factor(drop(y))]
  # x is required to be a matrix with more than more column
  x = as.matrix(x)
  if (!all(y %in% c(-1, 1))) 
    stop("y should be a factor with two levels.")
  nobs = as.integer(NROW(x))
  np = as.integer(NCOL(x))
  if (length(y) != nobs) 
    stop("x and y have different number of observations.")
  if (missing(kern)) {
    kern = rbfdot(sigma=sigest(x))
    cat("'kern' is missing: Gaussian kernel is used.\n")
  }
  ####################################################################
  # qval = as.double(qval)
  maxit = as.integer(maxit)
  eps = as.double(eps)
  gam = as.double(1e-7) # a tiny value to avoid matrix sigularity
  ####################################################################
  ## lambda setup
  if (missing(lambda)) {
      lambda = 10 ^ seq(0, -7, len=100)
    ulam = lambda
  } else {
    ulam = as.double(sort(lambda, decreasing=TRUE))
  }
  nlam = as.integer(length(lambda))
  fit = lhscpath(x, y, nobs, np, kern, ulam, nlam, eps, maxit, gam)
  fit.call = this.call
  class(fit) = c(class(fit), "lhsc")
  fit
} 

lhscpath = function(x, y, nobs, np, kern, ulam, nlam, eps, maxit, gam) {
  ####################################################################
    if (class(kern)[[1]] == "vanillakernel" && nobs >= np) {
    ## linear lhsc
      fit = .Fortran("llhsc", as.double(x),
        nobs, np, as.double(y), nlam, ulam, eps, maxit, gam,
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        btmat=double((np + 1) * nlam), PACKAGE="lhsc")
      anlam = fit$anlam
      alpha = matrix(fit$btmat[seq((np + 1) * anlam)], np + 1, anlam) 
    } else {
    ## kernel lhsc
      Kmat = kernelMatrix(kern, x)
      fit = .Fortran("klhsc", as.double(Kmat),
        nobs, as.double(y), nlam, ulam, eps, maxit, gam, 
        anlam=integer(1), npass=integer(nlam), jerr=integer(1), 
        alpmat=double((nobs + 1) * nlam), PACKAGE="lhsc")
      anlam = fit$anlam
      alpha = matrix(fit$alpmat[seq((nobs + 1) * anlam)], nobs + 1, anlam) 
    }

  ####################################################################
  ## wrap up output
  info = list(eps = eps, maxit = signif(maxit),
    kern = capture.output(show(kern)))
  outlist = list(alpha = alpha, lambda = ulam[seq(anlam)], 
    npass = fit$npass, jerr = fit$jerr, info = info)
  class(outlist) = c("lhscpath")
  outlist
}
