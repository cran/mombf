dmom <- function(x,V1=1,g=1,n=1,theta0,baseDensity='normal',nu=3) {

if (!(baseDensity %in% c('normal','t'))) stop("The only implemented baseDensity values are 'normal' and 't'")
if (missing(V1)) {
  if (is.vector(x)) V1 <- 1 else V1 <- diag(ncol(x))
}
if (missing(theta0)) {
  if (is.vector(V1)) theta0 <- 0 else theta0 <- rep(0,ncol(V1))
}
    
if (is.vector(V1)) {
  qtheta <- (x-theta0)^2/(n*g*V1)
  if (baseDensity=='normal') {
    ans <- qtheta*dnorm(x,theta0,sd=sqrt(n*g*V1))
  } else if (baseDensity=='t') {
    normct <- exp(lgamma(.5*(nu+1))-lgamma(.5*nu)-.5*log(nu*pi*n*g*V1))
    ans <- qtheta*normct*(1+qtheta/nu)^(-.5*(nu+1))*(nu-2)/nu
  }
} else {
  require(mvtnorm)
  qtheta <- mahalanobis(x,center=theta0,cov=n*g*V1)
  if (baseDensity=='normal') {
    ans <- qtheta*dmvnorm(x,mean=theta0,sigma=n*g*V1)/ncol(V1)
  } else if (baseDensity=='t') {
    ans <- qtheta*dmvt(x,delta=theta0,sigma=n*g*V1)*(nu-2)/(nu*length(theta0))
  }
}
return(ans)
}
