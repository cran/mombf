imomunknown <- function(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,ssr,B=10^5) {
#Inverse moment Bayes factor for linear models (unknown sigma^2 case). Only implements case for prior parameter k=1.
# - theta1hat: vector with estimated value of the coefficients that are to be tested
# - V1: submatrix of covariance corresponding to elements in theta1hat
# - n: sample size used to fit the model
# - nuisance.theta: number of nuisance regression coefficients (including intercept)
# - g: prior parameter
# - nu: prior parameter governing the tail behavior. nu=1 specifies Cauchy-like tails.
# - theta0: hypothesized value for theta1hat (defaults to 0)
# - ssr: sum of squared residuals
# - B: number of Monte Carlo samples to approximate the integral
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat)
sigma2 <- 1/rgamma(B,(n-nuisance.theta)/2,ssr/2)
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1) / sigma2 #noncentr param
z <- rchisq(B,df=p1,ncp=l)
m <- mean((n*g/z)^((nu+p1)/2) * exp(-n*g/z) )
t1 <- theta1hat-theta0; t1 <- matrix(t1,nrow=1) %*% solve(V1) %*% matrix(t1,ncol=1) / ssr #noncentr param
bf <- exp((p1/2)*log(2/(n*g)) + lgamma(p1/2) - lgamma(nu/2) + ((n-nuisance.theta)/2)*log(1+t1)) * m
return(bf)
}
