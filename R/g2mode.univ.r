g2mode.univ <- function(g,V1,n,prior='iMom') {
#Find Mom/iMom prior mode for |(theta-theta0)/sigma| for a given vector of the prior parameter g values
  if (prior=='Mom') {
    return(sqrt(2*n*V1*g))
  } else if (prior=='iMom') {
    return(sqrt(n*V1*g))
  }
}
