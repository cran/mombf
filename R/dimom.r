dimom <- function(x,theta0,V1,g=1,n=1,nu=1,logscale=FALSE) {
# x: vector or matrix with values where the density should be evaluated (each row in the matrix corresponds to a new observation)
# theta0: location parameter
# V1: prior covariance parameter (matrix)
# g: prior variance parameter (scalar)
# n: prior variance parameter (scalar). Actually only the value of n*g is taken into account
if (is.vector(x)) {
  if (missing(theta0)) theta0 <- 0
  if (missing(V1)) V1 <- 1
  qtheta <- (x-theta0)^2/(n*g*V1)
  k <- -0.5*log(sqrt(n*g*V1)) - lgamma(nu/2)
  p1 <- 1
} else {
  if (missing(theta0)) theta0 <- rep(0,ncol(x))
  if (missing(V1)) V1 <- diag(ncol(x))
  qtheta <- t(matrix(x,nrow=nrow(x))) - theta0
  qtheta <- qtheta %*% solve(n*g*V1) %*% t(qtheta)
  k <- -0.5*log(sqrt(n*g*det(V1))) - lgamma(nu/2) + lgamma(ncol(x)/2) - .5*ncol(x)*log(pi)
  p1 <- ncol(x)
}
ans <- (k - .5*(nu+p1)*log(qtheta) -1/qtheta)
if (!logscale) { ans <- exp(ans); ans[is.na(ans)] <- 0 }
return(ans)
}
