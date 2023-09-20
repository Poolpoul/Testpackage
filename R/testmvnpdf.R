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



