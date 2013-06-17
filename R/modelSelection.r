###
### modelSelection.R
###

### Methods for msfit objects

setMethod("show", signature(object='msfit'), function(object) {
  cat('msfit object with',ncol(object$postSample),'variables\n')
  ifelse(any(object$postMode!=0), paste('  Posterior mode: covariate',which(object$postMode==1)), '  Posterior mode: null model')
  cat("Use postprob() to get posterior model probabilities\n")
  cat("Elements $margpp, $postMode, $postSample and $coef contain further information (see help('msfit') and help('modelSelection') for details)\n")
}
)

setMethod("postProb", signature(object='msfit'), function(object, nmax) {
  modelpp <- unique(data.frame(object$postSample==1, logpp=object$postProb))
  modelpp <- data.frame(modelid= apply(modelpp[,1:(ncol(modelpp)-1)], 1, function(z) paste(which(z),collapse=',')), logpp=modelpp$logpp)
  modelpp$logpp <- modelpp$logpp - modelpp$logpp[1]
  modelpp$pp <- exp(modelpp$logpp)/sum(exp(modelpp$logpp))
  modelpp <- modelpp[order(modelpp$pp,decreasing=TRUE),]
  if (!missing(nmax)) modelpp <- modelpp[1:nmax,]
  modelpp[,c('modelid','pp')]
}
)



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
  } else if (method=='auto') {
    if (priorCoef@priorDistr!='pMOM') { 
      warning("method=='auto' is only available for 'pMOM' priors. Using method=='Laplace' instead")
      method <- as.integer(-1)
    } else {
      method <- as.integer(2)
    }
  } else {
    stop('Invalid argument method')
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
    parprDeltap <- double(2)
  } else if (priorDelta@priorDistr=='binomial') {
    if ('p' %in% priorDelta@priorPars) {
      prDelta <- as.integer(1)
      prDeltap <- as.double(priorDelta@priorPars['p'])
      if ((prDeltap<=0) | (prDeltap>=1)) stop("p must be between 0 and 1 for priorDelta@priorDistr=='binomial'")
      parprDeltap <- double(2)
    } else {
      prDelta <- as.integer(2)
      prDeltap <- as.double(.5)
      parprDeltap <- as.double(priorDelta@priorPars[c('alpha.p','beta.p')])
    }
  } else {
    stop('Prior specified in priorDelta not recognized')
  }

  #Initialize
  postMode <- rep(as.integer(0),p); postModeProb <- double(1)
  if (initSearch=='greedy') {
    niterGreed <- as.integer(100)
    ans <- .Call("greedyVarSelCI", postMode,postModeProb,knownphi,prior,niterGreed,ndeltaini,deltaini,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,phi,tau,r,prDelta,prDeltap,parprDeltap,as.integer(verbose))
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
  if (prDelta==2) postOther <- double(mcmc2save) else postOther <- double(0)
  margpp <- double(p); postProb <- double(mcmc2save)
  ans <- .Call("modelSelectionCI", postSample,postOther,margpp,postMode,postModeProb,postProb,knownphi,prior,niter,thinning,burnin,ndeltaini,deltaini,n,p,y,sumy2,as.double(x),XtX,ytX,method,B,alpha,lambda,phi,tau,r,prDelta,prDeltap,parprDeltap,as.integer(verbose))
  postSample <- matrix(postSample,ncol=p)
  coef <- rep(0,ncol(x))
  if (any(postMode)) {
    pm <- postMode(y=y,x=x[,postMode==1,drop=FALSE],priorCoef=priorCoef)
    coef[postMode==1] <- pm$coef
  }
  ans <- list(postSample=postSample,postOther=postOther,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb,coef=coef)
  new("msfit",ans)
}

modelSelectionR <- function(y, x, niter=10^4, marginalFunction, priorFunction, betaBinPrior, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - niter: number of Gibbs sampling iterations
# - marginalFunction: function to compute the marginal density of the data under each model
# - priorFunction: function to compute the model prior probability
# - betaBinPrior: if specified, priorFunction argument is ignored and set to a binomial prior with Beta-hyperprior for the success prob. betaBinPrior should be a vector with named elements 'alpha.p' and 'beta.p', e.g. betaBinPrior= c(alpha.p=1,beta.p=1)
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# ...: other arguments to be passed on to marginalFunction
# Output: list
# - postSample: posterior samples for model indicator
# - postOther: posterior samples for other parameters (probBin: success probability for Binomial prior on number of coefficients in the model)
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)
p <- ncol(x)
if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
if (!missing(betaBinPrior)) {
  #Initialize probBin
  if ((betaBinPrior['alpha.p']>1) && (betaBinPrior['beta.p']>1)) {
    probBin <- (betaBinPrior['alpha.p']-1)/(betaBinPrior['alpha.p']+betaBinPrior['beta.p']-2)
  } else {
    probBin <- (betaBinPrior['alpha.p'])/(betaBinPrior['alpha.p']+betaBinPrior['beta.p'])
  }
  postOther <- matrix(NA,nrow=niter,ncol=1); colnames(postOther) <- c('probBin')
  priorFunction <- function(sel, logscale=TRUE) dbinom(x=sum(sel),size=length(sel),prob=probBin,log=logscale)
} else {
  postOther <- matrix(NA,nrow=niter,ncol=0)
}
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
  if (!missing(betaBinPrior)) {
    probBin <- rbeta(1,betaBinPrior['alpha.p']+sum(sel), betaBinPrior['beta.p']+sum(!sel))
    postOther[i,'probBin'] <- probBin
  }
  postSample[i,] <- sel
  postProb[i] <- currentM + currentP
  if (verbose && ((i%%niter10)==0)) cat('.')
}
if (verbose) cat('\n')
ans <- list(postSample=postSample,postOther=postOther,margpp=margpp/niter,postMode=postMode,postModeProb=postModeProb,postProb=postProb)
#new("msfit",ans)
return(ans)
}


#Gibbs model selection using BIC to approximate marginal likelihood
# - y, x, xadj: response, covariates under selection and adjustment covariates
# - family: glm family, passed on to glm
# - niter: number of Gibbs iteration
# - burnin: burn-in iterations
# - modelPrior: function evaluating model log-prior probability. Takes a logical vector as input
# Returns:
# - postModel, postCoef1, postCoef2, margpp (analogous to pmomPM. postCoef1 & postCoef2 are MLEs under each visited model)
modelselBIC <- function(y, x, xadj, family, niter=1000, burnin= round(.1*niter), modelPrior, verbose=TRUE) {
  pluginJoint <- function(sel) {
    p <- sum(sel)
    ans <- vector("list",2); names(ans) <- c('marginal','coef')
    if (p>0 & p<=length(y)) {
      glm1 <- glm(y ~ x[,sel,drop=FALSE] + xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- coef(glm1)[1:p]; ans$coef2 <- coef(glm1)[-1:-p]
    } else if (p==0) {
      glm1 <- glm(y ~ xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- double(0); ans$coef2 <- coef(glm1)
    } else { ans$marginal <- -Inf; ans$coef1 <- ans$coef2 <- 0}
    return(ans)
  }
  #Greedy iterations
  if (verbose) cat("Initializing...")
  sel <- rep(FALSE,ncol(x))
  mcur <- pluginJoint(sel)$marginal
  nchanges <- 1; it <- 1
  while (nchanges>0 & it<100) {
    nchanges <- 0; it <- it+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- pluginJoint(selnew)$marginal
      if (mnew>mcur) { sel[i] <- selnew[i]; mcur <- mnew; nchanges <- nchanges+1 }
    }
  }
  if (verbose) { cat(" Done\nGibbs sampling") }
  #Gibbs iterations
  niter10 <- ceiling(niter/10)
  postModel <- matrix(NA,nrow=niter,ncol=ncol(x))
  postCoef1 <- matrix(0,nrow=niter,ncol=ncol(x))
  postCoef2 <- matrix(0,nrow=niter,ncol=ncol(xadj))
  curmod <- pluginJoint(sel)
  for (j in 1:niter) {
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      newmod <- pluginJoint(selnew)
      mnew <- newmod$marginal
      if (runif(1) < 1/(1+exp(mcur-mnew))) { sel[i] <- selnew[i]; mcur <- mnew; curmod <- newmod }
    }
    postModel[j,] <- sel
    postCoef1[j,sel] <- curmod$coef1; postCoef2[j,] <- curmod$coef2
    if (verbose & ((j%%niter10)==0)) cat(".")
  }
  if (verbose) cat("Done\n")
  #Return output
  if (burnin>0) { postModel <- postModel[-1:-burnin,,drop=FALSE]; postCoef1 <- postCoef1[-1:-burnin,,drop=FALSE]; postCoef2 <- postCoef2[-1:-burnin,,drop=FALSE] }
  ans <- list(postModel=postModel, postCoef1=postCoef1, postCoef2=postCoef2, margpp=colMeans(postModel))
}


## Common prior distributions on model space

binomPrior <- function(sel, prob=.5, logscale=TRUE) {  dbinom(x=sum(sel),size=length(sel),prob=prob,log=logscale) }
unifPrior <- function(sel, logscale=TRUE) { ifelse(logscale,-length(sel)*log(2),2^(-length(sel)))  }
bbPrior <- function(sel, alpha=1, beta=1, logscale=TRUE) {
  ans <- lbeta(sum(sel) + alpha, sum(!sel) + beta) - lbeta(alpha,beta)
  ifelse(logscale,ans,exp(ans))
}
