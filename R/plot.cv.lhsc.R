plot.cv.lhsc = function(x, sign.lambda=1, ...) {
  cvobj = x
  xlab = "log(Lambda)"
  if (sign.lambda < 0) 
    xlab = paste0("-", xlab)
  plot.args = list(x = sign.lambda * log(cvobj$lambda), 
    y = cvobj$cvm, ylim = range(cvobj$cvup, cvobj$cvlo), 
    xlab = xlab, ylab = cvobj$name, type = "n")
  new.args = list(...)
  if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, 
    cvobj$cvlo, width=0.01, col="darkgrey")
  points(sign.lambda * log(cvobj$lambda), 
    cvobj$cvm, pch=20, col="red")
  axis(side=3, at=sign.lambda * log(cvobj$lambda), 
    labels=paste(cvobj$nz), tick=FALSE, line=0)
  abline(v=sign.lambda * log(cvobj$lambda.min), lty=3)
  abline(v=sign.lambda * log(cvobj$lambda.1se), lty=3)
  invisible()
} 
