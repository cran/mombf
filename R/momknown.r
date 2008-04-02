momknown <- function(theta1hat,V1,n,g=1,theta0,sigma,logbf=FALSE) {
#Moment Bayes factor for linear models (known sigma^2 case). Only implements case for prior parameter k=1.
# - theta1hat: vector with estimated value of the coefficients that are to be tested
# - V1: submatrix of covariance corresponding to elements in theta1hat
# - n: sample size used to fit the model
# - g: prior parameter
# - theta0: hypothesized value for theta1hat (defaults to 0)
# - sigma: residual standard deviation, assumed to be known
# - logbf: if log==TRUE, the natural logarithm of the Bayes factor is returned
if (missing(sigma)) stop('sigma must be specified')
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat)
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) * n*g/((1+n*g)*sigma^2) #noncentr param
muk <- p1+l
t1 <- matrix(theta1hat-theta0,nrow=1) %*% solve(V1) %*% matrix(theta1hat-theta0,ncol=1) * n*g/((1+n*g)*sigma^2)
bf <- .5*t1 + log(muk) - log(1+n*g) - (p1/2)*log(1+n*g)
if (!logbf) bf <- exp(bf)
return(bf)
}
