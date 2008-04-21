pmom <- function(q,V1=1,g=1,n=1) {

  z <- .5-(pnorm(abs(q)/sqrt(n*V1*g)) - abs(q)/sqrt(2*pi*n*V1*g) * exp(-.5*q^2/(g*n*V1)) - .5)
  return(z*(q<=0)+(1-z)*(q>0))
}
