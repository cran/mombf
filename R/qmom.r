qmom <- function(p,V1=1,g=1,n=1) {

  e <- function(q) { return((pmom(q,V1,g,n)-pneg)^2) }
  ans <- double(length(p))
  for (i in 1:length(p)) {
    pneg <- ifelse(p[i]<=.5,p[i],1-p[i])
    ans[i] <- nlminb(start=-1,objective=e)$par
    if (p[i]>.5) ans[i] <- -ans[i]
  }
  return(ans)
}
