dmom <- function(x,V1=1,g=1,n=1,theta0) {

if (missing(V1)) {
  if (is.vector(x)) V1 <- 1 else V1 <- diag(ncol(x))
}
if (missing(theta0)) {
  if (is.vector(V1)) theta0 <- 0 else theta0 <- rep(0,ncol(V1))
}
    
if (is.vector(V1)) {
  qtheta <- (x-theta0)^2/(n*g*V1)
  return(qtheta*dnorm(x,theta0,sd=sqrt(n*g*V1)))
} else {
  require(mvtnorm)
  qtheta <- mahalanobis(x,center=theta0,cov=n*g*V1)
  return(qtheta*dmvnorm(x,mean=theta0,sigma=n*V1)/ncol(V1))
}
}
