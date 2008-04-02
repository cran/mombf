g2mode <- function(g,n,prior='iMom',nu=1,p) {
#Find Mom/iMom prior mode for (theta-theta0)'(X'X)(theta-theta0)/(sigma*n) for a given vector of the prior parameter g values
  if ((missing(p)) & (prior='iMom')) stop('p must be specified for prior=="iMom"')
  if (p==1) cat('For p==1 (univariate setting) it may be more intuitive to use the function mode2g.univ\n')
  if (prior=='Mom') {
    return(2*g)
  } else if (prior=='iMom') {
    return(sqrt(2*g/(nu+p)))
  }
}
