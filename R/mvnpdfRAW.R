#' Title
#'
#' description
#'
#' details
#'
#'
#' @param x x is a matrix containing the observation in columns
#' @param mean mean vector of the multivariate Gaussian distribution
#' @param varcovM variance covariance matrix
#' @param Log logical parameter (default = TRUE)
#'
#' @return value -- list of the results
#' @export
#'
#' @examples
#' X <- matrix(c(-0.5,1.5,0,1, -1, 1), nrow = 2)
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

#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from a call of
#' \code{mnvpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} function.
#'
#' @return Nothing is returned, only a plot is given.
#' @export
#'
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by = 0.1), nrow = 1), Log=FALSE)
#' plot(pdfvalues)
plot.mvnpdf <- function(x, ...) {
  plot(x$x, x$y, type = "l", ...)
}

