
#### General model selection routines
modelSelection <- function(y, x, center=TRUE, scale=TRUE, niter=10^4, thinning=1, burnin=round(niter/10), priorCoef, priorDelta, priorVar, phi, deltaini=rep(FALSE,ncol(x)), initSearch='greedy', method, B=10^5, verbose=TRUE) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - center: if center==TRUE y and x are centered to have zero mean, therefore eliminating the need to include an intercept term in x.
# - scale: if scale==TRUE y and columns in x are scaled to have standard deviation 1
# - niter: number of Gibbs sampling iterations
# - thinning: MCMC thinning factor, i.e. only one out of each thinning iterations are reported. Defaults to thinning=1, i.e. no thinning
# - burnin: number of burn-in MCMC iterations. Defaults to 10% of niter. Set to 0 for no burn-in.
# - priorCoef: prior distribution for the coefficients. Must be object of class 'msPriorSpec' with slot priorType set to 'coefficients'. Possible values for slot priorDistr are 'pMOM', 'piMOM' and 'peMOM'.
# - priorDelta: prior on model indicator space. Must be object of class 'msPriorSpec' with slot priorType set to 'modelIndicator'. Possible values for slot priorDistr are 'uniform' and 'binomial'
# - priorVar: prior on residual variance. Must be object of class 'msPriorSpec' with slot priorType set to 'nuisancePars'. Slot priorDistr must be equal to 'invgamma'.
# - phi: residual variance. Typically this is unknown and therefore left missing. If specified argument priorVar is ignored.
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# - initSearch: algorithm to refine deltaini. initSearch=='greedy' uses a greedy Gibbs sampling search. initSearch=='SCAD' sets deltaini to the non-zero elements in a SCAD fit with cross-validated regularization parameter. initSearch=='none' leaves deltaini unmodified.
# - method: method to compute marginal densities. method=='Laplace' for Laplace approx, method=='MC' for Importance Sampling, method=='Hybrid' for Hybrid Laplace-IS (the latter method is only used for piMOM prior with unknown residual variance phi)
# - B: number of samples to use in Importance Sampling scheme. Ignored if method=='Laplace'.
# - verbose: set verbose==TRUE to print iteration progress
# Output: list
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)

  #Check input
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  ct <- (colMeans(x^2)-colMeans(x)^2)==0
  y <- scale(y,center=center,scale=scale); x[,!ct] <- scale(x[,!ct],center=center,scale=scale)
  if (missing(priorCoef)) { priorCoef <- new("msPriorSpec",priorType='coefficients',priorDistr='pMOM',priorPars=c(tau=1,r=1)) }
  if (missing(priorDelta)) { priorDelta <- new("msPriorSpec",priorType='modelIndicator',priorDistr='uniform',priorPars=double(0)) }
  if (missing(priorVar)) { priorVar <- new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=.01,lambda=.01)) }
  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) } 
  p <- ncol(x); n <- length(y)
  if (missing(deltaini)) { 
    deltaini <- integer(0); ndeltaini= as.integer(0) 
  } else { 
    if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
    if (!is.logical(deltaini)) { stop('deltaini must be of type logical') } else { ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(which(deltaini)-1) }
  }
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (method=='Laplace') {
    method <- as.integer(0)
  } else if (method=='MC') {
    method <- as.integer(1)
  } else if (method=='Hybrid') {
    if ((priorCoef@priorDistr!='piMOM') | (knownphi==1)) { 
      warning("method=='Hybrid' is only available for 'piMOM' priors with unknown phi. Using method=='Laplace' instead")
      method <- as.integer(0)
    } else {
      method <- as.integer(2)
    }
  }

  #Format arguments for .Call
  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning); B <- as.integer(B)
  sumy2 <- as.double(sum(y^2)); XtX <- t(x) %*% x; ytX <- as.vector(matrix(y,nrow=1) %*% x)

  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r']); prior <- as.integer(0)
  } else if (priorCoef@priorDistr=='piMOM') {
    r <- as.integer(1); prior <- as.integer(1)
  } else if (priorCoef@priorDistr=='peMOM') {
    r <- as.integer(1); prior <- as.integer(2)
  } else {
    stop('Prior specified in priorDistr not recognized')
  }
  tau <- as.double(priorCoef@priorPars['tau'])
  alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda']) 
  if (priorDelta@priorDistr=='uniform') {
    prDelta <- as.integer(0)
    prDeltap <- as.double(0)
  } else if (priorDelta@priorDist=='binomial') {
    prDelta <- as.integer(1)
    prDeltap <- as.double(priorDelta@priorPars['p'])
    if ((prDeltap<=0) | (prDeltap>=1)) stop("p must be between 0 and 1 for priorDelta@priorDist=='binomial'")
  } else {
    stop('Prior specified in priorDelta not recognized')
  }

  #Initialize
  postMode <- rep(as.integer(0),p); postModeProb <- double(1)
  if (initSearch=='greedy') {
    niterGreed <- as.integer(100)
    ans <- .Call("greedyVarSelCI", postMode,postModeProb,knownphi,prior,niterGreed,ndeltaini,deltaini,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,phi,tau,r,prDelta,prDeltap,as.integer(verbose))
    ndeltaini <- as.integer(sum(postMode)); deltaini <- as.integer(which(as.logical(postMode))-1)
  } else if (initSearch=='SCAD') {
    require(ncvreg)
    if (verbose) cat("Initializing via SCAD cross-validation...")
    deltaini <- rep(TRUE,ncol(x))
    cvscad <- cv.ncvreg(X=x[,!ct],y=y-mean(y),family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
    deltaini[!ct] <- ncvreg(X=x[,!ct],y=y-mean(y),penalty='SCAD',dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1,1]!=0
    ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(which(deltaini)-1)
    if (verbose) cat(" Done\n")
  }
  
  #Run MCMC
  mcmc2save <- floor((niter-burnin)/thinning)
  postSample <- rep(as.integer(0),p*mcmc2save)
  margpp <- double(p); postProb <- double(mcmc2save)
  ans <- .Call("modelSelectionCI", postSample,margpp,postMode,postModeProb,postProb,knownphi,prior,niter,thinning,burnin,ndeltaini,deltaini,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,phi,tau,r,prDelta,prDeltap,as.integer(verbose))
  postSample <- matrix(postSample,ncol=p)
  return(list(postSample=postSample,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb))
}

modelSelectionR <- function(y, x, niter=10^4, marginalFunction, priorFunction, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - niter: number of Gibbs sampling iterations
# - marginalFunction: function to compute the marginal density of the data under each model
# - priorFunction: function to compute the model prior probability
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# ...: other arguments to be passed on to marginalFunction
# Output: list
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)
p <- ncol(x)
if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
sel <- postMode <- deltaini
currentM <- marginalFunction(y=y,x=x[,sel,drop=FALSE],logscale=TRUE,...)
currentP <- priorFunction(sel,logscale=TRUE)
postModeProb <- currentM + currentP
postSample <- matrix(NA,nrow=niter,ncol=p)
margpp <- double(p)
postProb <- double(niter)
postProb[1] <- postModeProb
niter10 <- ceiling(niter/10)
for (i in 1:niter) {
  for (j in 1:p) {
    selnew <- sel; selnew[j] <- !sel[j]
    newM <- marginalFunction(y=y,x=x[,selnew,drop=FALSE],logscale=TRUE,...)
    newP <- priorFunction(selnew,logscale=TRUE)
    if (newM+newP>postModeProb) {
      postModeProb <- newM+newP
      postMode <- selnew
    }
    pp <- 1/(1+exp(-currentM+newM-currentP+newP))
    if (sel[j]) {  #if variable in the model
      sel[j] <- runif(1)<pp
      margpp[j] <- margpp[j]+pp
    } else {       #if variable not in the model
      sel[j] <- runif(1)>pp
      margpp[j] <- margpp[j]+1-pp
    }
    if (sel[j]==selnew[j]) {  #if value was updated, update marginal and prior densities
      currentM <- newM
      currentP <- newP
    }
  }
  postSample[i,] <- sel
  postProb[i] <- currentM + currentP
  if (verbose && ((i%%niter10)==0)) cat('.')
}
if (verbose) cat('\n')
return(list(postSample=postSample,margpp=margpp/niter,postMode=postMode,postModeProb=postModeProb,postProb=postProb))
}
