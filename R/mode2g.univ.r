mode2g.univ <- function(prior.mode,V1,n,prior='iMom') {
#Find the value of the prior parameter g value giving a certain Mom or iMom prior mode for |(theta-theta0)/sigma|
  if (prior=='Mom') {
    return(prior.mode^2/(2*n*V1))
  } else if (prior=='iMom') {
    return(prior.mode^2/(n*V1))
  }
}
