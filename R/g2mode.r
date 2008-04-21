g2mode <- function(g,prior='iMom',nu=1,dim=1) {
  if (prior=='Mom') {
    return(2*g)
  } else if (prior=='iMom') {
    return(2*g/(nu+dim))
  }
}
