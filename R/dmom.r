#Wrapper to call dpmom (product mom) and dqmom (quadratic mom)
dmom <- function(x, tau=1, phi=1, r=1, V1, baseDensity='normal', nu=3, logscale=FALSE, penalty='product') {
  if (penalty=='product') {
    ans <- dpmom(x, tau=tau, phi=phi, r=r, baseDensity=baseDensity, logscale=logscale)
  } else if (penalty=='quadratic') {
    ans <- dqmom(x, V1=V1, g=tau, n=1, baseDensity=baseDensity, nu=nu)
  } else {
    stop("Only 'penalty==product' and 'penalty==quadratic' are implemented")
  }
  return(ans)
}

##Product MOM density
setGeneric("dpmom", function(x, tau=1, phi=1, r=1, baseDensity='normal', logscale=FALSE) standardGeneric("dpmom"))
setMethod("dpmom", signature(x='vector'), function(x, tau=1, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
  if (baseDensity!='normal') stop("Only baseDensity=='normal' is implemented for the product MOM")
  ans <- dnorm(x,0,sd=sqrt(tau*phi),log=TRUE) + r*log(x^2/(tau*phi)) - sum(log(seq(1,2*r-1,by=2)))
  if (!logscale) ans <- exp(ans)
  ans
}
)
setMethod("dpmom", signature(x='matrix'), function(x, tau=1, phi=1, r=1, baseDensity='normal', logscale=FALSE) {
  if (baseDensity!='normal') stop("Only baseDensity=='normal' is implemented for the product MOM")
  p <- ncol(x)
  normct <- p*sum(log(seq(1,2*r-1,by=2)))
  distval <- rowSums(x^2)
  ans <- -(p * log(2 * pi) + p*(log(phi)+log(tau))  + distval/(phi*tau))/2 + r*rowSums(log(x^2/(tau*phi))) - normct
  if (!logscale) ans <- exp(ans)
  ans
}
)

##Quadratic MOM density
dqmom <- function(x,V1=1,g=1,n=1,theta0,baseDensity='normal',nu=3) {

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
