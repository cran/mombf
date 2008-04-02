mombf.lm <- function(lm1,coef,g,prior.mode,theta0,logbf=FALSE) {
#Moment Bayes factor for linear models. Wrapper to make it easier to call momunknown
#Input:
# - lm1: object of type lm
# - coef: vector with indexes of coefficients to be tested. e.g. coef=c(2,3) and theta0=c(0,0) tests coef(lm1)[2]=coef(lm1)[3]=0
# - g: vector with prior parameter values
# - prior.mode: if specified, g is determined by calling g2mode. prior.mode is the mode a priori for (theta-theta0)'(X'X)(theta-theta0)/(n*sigma^2)
# - theta0: hypothesized value for coef(lm1)[coef] (defaults to 0)
#Output: Bayes factor based on Moment prior
  
if ((!missing(g)) & (!missing(prior.mode))) warning('Both g and prior.mode were specified. g will be ignored')
if ((missing(g)) & (missing(prior.mode))) stop('Either g or prior.mode must be specified')
if (missing(theta0)) theta0 <- rep(0,length(coef)) else if (length(theta0)!=length(coef)) stop('theta0 must have the same length as coef')
  
  thetahat <- coef(lm1)
  V <- summary(lm1)$cov.unscaled
  n <- length(lm1$residuals); p <- length(thetahat); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) stop('Non-valid value for coef. Use only values between 1 and the number of coefficients in lm1')
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  if (!missing(prior.mode)) g <- mode2g(prior.mode,prior='Mom')
  bf.mom <- momunknown(thetahat[coef],V[coef,coef],n=n,g=g,ssr=ssr,nuisance.theta=p-p1,theta0=theta0,logbf=logbf)
  return(bf.mom)
}
