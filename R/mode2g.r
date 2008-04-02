mode2g <- function(prior.mode,prior='iMom',nu=1,p) {
#Find the prior parameter g value given either the Mom or iMom prior mode for (theta-theta0)'(X'X)(theta-theta0)/(sigma*n)
  if ((missing(p)) & (prior='iMom')) stop('p must be specified for prior=="iMom"')
  if (p==1) cat('For p==1 (univariate setting) it may be more intuitive to use the function mode2g.univ\n')
  if (prior=='Mom') {
    return(prior.mode^2/2)
  } else if (prior=='iMom') {
    return(prior.mode^2 *(nu+p)/2)
  }
}
