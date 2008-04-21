dmom <- function(x,V1=1,g=1,n=1,theta0) {
  
if (is.vector(x)) {
  if (missing(theta0)) theta0 <- 0
  if (missing(V1)) V1 <- 1
  qtheta <- (x-theta0)^2/(n*g*V1)
  return(qtheta*dnorm(x,theta0,sd=sqrt(n*g*V1)))
} else {
  require(mvtnorm)
  if (missing(theta0)) theta0 <- rep(0,ncol(x))
  if (missing(V1)) V1 <- diag(ncol(x))
  qtheta <- t(matrix(x,nrow=nrow(x))) - theta0
  qtheta <- qtheta %*% solve(n*g*V1) %*% t(qtheta)
  return(qtheta*dmvnorm(x,mean=theta0,sigma=V1)/ncol(x))
}
}
