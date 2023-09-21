#microbenchmark
n <- 1000
mb <- microbenchmark(mvtnorm::dmvnorm(matrix(1.96, nrow = n, ncol = 2)),
                     mvnpdf(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE),
                     mvnpdfsmart(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE),
                     mvnpdfoptim(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE),
                     times=100L)
mb

#mvnpdf
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))
  
  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }
  
  if (!Log) {
    y <- exp(y)
  }
  
  res <- list(x = x, y = y)
  return(res)
}

#mvnpdfsmart
mvnpdfsmart <- function(x, mean =  rep(0, nrow(x)),
                        varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))
  
  y <- rep(NA, n)
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y[j] <- yj
  }
  
  if (!Log) {
    y <- exp(y)
  }
  
  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}

##opti
mvnpdfoptim <- function(x, mean =  rep(0, nrow(x)),
                        varcovM = diag(nrow(x)), Log=TRUE){
  
  if(!is.matrix(x)){
    x <- matrix(x, ncol=1)
  }
  
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x-mean
  
  Rinv <- backsolve(chol(varcovM), x=diag(p))
  xRinv <- apply(X=x0, MARGIN=2, FUN=crossprod, y=Rinv)
  logSqrtDetvarcovM <- sum(log(diag(Rinv)))
  
  quadform <- apply(X=xRinv, MARGIN=2, FUN=crossprod)
  y <- (-0.5*quadform + logSqrtDetvarcovM -p*log(2*pi)/2)
  
  if(!Log){
    y <- exp(y)
  }
  
  res <- list(x = x, y = y)
  class(res) <- "mvnpdf"
  return(res)
}

#smart -> select lines -> Profile -> profile selected line -> stop profiling when finished
n <- 10e4
pdfval <- mvnpdfsmart(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE)

#mvnpdf
n <- 10e4
pdfval <- mvnpdf(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE)

#opti
n <- 10e4
profvis::profvis(mvnpdfoptim(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE))

mvnpdfoptim(x=matrix(1.96, nrow = 2, ncol = n), Log=FALSE)


#################################################################################
#RCPP
