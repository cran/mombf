mode2g <- function(prior.mode,prior='iMom',nu=1,dim=1) {
  if (prior=='Mom') {
    return(prior.mode/2)
  } else if (prior=='iMom') {
    return(prior.mode *(nu+dim)/2)
  }
}
