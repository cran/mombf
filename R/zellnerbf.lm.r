zellnerbf.lm <- function(lm1,coef,g,theta0,logbf=FALSE) {
#Bayes factor based on Zellner's g-prior for linear models. Wrapper to make it easier to call zbfunknown
#Input:
# - lm1: object of type lm
# - coef: vector with indexes of coefficients to be tested. e.g. coef=c(2,3) and theta0=c(0,0) tests coef(lm1)[2]=coef(lm1)[3]=0
# - g: prior for coef(lm1)[coef] under the alternative is N(theta0,g*sigma^2*n*(X'X))
# - theta0: hypothesized value for coef(lm1)[coef] (defaults to 0)
#Output: Bayes factor based on Moment prior
  
if (missing(g)) stop('g must be specified')
if (missing(theta0)) theta0 <- rep(0,length(coef)) else if (length(theta0)!=length(coef)) stop('theta0 must have the same length as coef')
  
  thetahat <- coef(lm1)
  V <- summary(lm1)$cov.unscaled
  n <- length(lm1$residuals); p <- length(thetahat); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) stop('Non-valid value for coef. Use only values between 1 and the number of coefficients in lm1')
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  bf.zellner <- zbfunknown(thetahat[coef],V[coef,coef],n=n,nuisance.theta=p-p1,g=g,theta0=theta0,ssr=ssr,logbf=logbf)
  return(bf.zellner)
}
