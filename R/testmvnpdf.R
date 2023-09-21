#Créer un matrice
M <- matrix(c(1,2,2,3), nrow = 2, byrow = TRUE)

#det = determinant
det(M)
x <- matrix(c(-0.5, 1.5), nrow = 2)
x
#solv M = inverse de la matrix de M
#Multiplier des matrices entre elles
solve(M)
M%*% solve(M)

t(x)
t(x) %*% solve(M) %*% x
# bouclefort
f <- function(x) {
  return(x^2)
}

##
y <- runif(10)
#c = concatener
#z = result
z <- c()
for (i in 1:10) {
 z <- c(z, f(y[i]))
}
z

#mvnpdf the density function for a multivariateGaussian distribution is:
mvnpdftest <- function(x, mean, varcovM, Log =TRUE){
#nombre de ligne de x -> P
p = nrow(varcovM)
(2*pi)^(-p/2)*det(varcovM)^(1/2)*exp((-1/2)*t(x-mean)%*%solve(varcovM)%*%(x-mean))
}
#transposer matrix : t()

x <- matrix(c(-0.5,1.5), nrow = 2)
mean <- matrix(c(0,0), nrow = 2)
varcovM <- matrix(c(1,0,0,1), nrow = 2)
X <- matrix(c(-0.5,1.5,0,1, -1, 1), nrow = 2)
mvnpdftest(x, mean, varcovM)
mvtnorm::dmvnorm(t(x), mean, varcovM)

#ajouter une boucle fort (calcul pour chaque column)
results <- c()
for (j in 1:ncol(X)) {
  results <- c(results, mvnpdftest(X[,j,drop = FALSE],mean, varcovM))
}
results
#comparer mvnpdf et mvnpdfsmart
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
