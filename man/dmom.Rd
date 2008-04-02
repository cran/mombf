\name{dmom}
\alias{dmom}
\alias{dimom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Density for Moment prior and inverse moment prior. }
\description{
  These functions evaluate the density for the moment and inverse moment
  priors, in the univariate and multivariate setting.
}
\usage{
dmom(x, theta0, V1, g = 1, n = 1)
dimom(x, theta0, V1, g = 1, n = 1, nu = 1, logscale = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ In the univariate setting, \code{x} is a vector with the
    values at which to evaluate the density. In the multivariate setting
  it is a matrix with an observation in each row.}
  \item{theta0}{ Location parameter. Defaults to 0. }
  \item{V1}{ Scale matrix. Defaults to 1 in univariate setting and
    the identity matrix in the multivariate setting. }
  \item{g}{ Prior parameter. See details. }
  \item{n}{ Prior parameter. See details. }
  \item{nu}{ Prior parameter. The
    tails of the inverse moment prior are proportional to the tails of a
    multivariate T with \code{nu} degrees of freedom.}
  \item{logscale}{ For \code{logscale==TRUE}, \code{dimom} returns the
    natural log of the prior density.} 
}
\details{
 Define the quadratic form \code{q(theta)= (theta-theta0)' *
   solve(V1) * (theta-theta0) / (n*g)}.
   The moment prior density is proportional to
   q(theta)*dmvnorm(theta,theta0,n*g*V1).
   The inverse moment prior density is proportional to
   \code{q(theta)^(-(nu+d)/2) * exp(-1/q(theta))}.
}
\value{
  \code{dmom} returns the value of the moment prior density.
  \code{dimom} returns the value of the inverse moment prior density.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\seealso{ \code{\link{g2mode}} and \code{\link{g2mode.univ}} to find the
  prior mode corresponding to a given \code{g}. \code{\link{mode2g}} and
  \code{\link{mode2g.univ}} to find the g value corresponding to a given
prior mode.}
\examples{
#evaluate and plot the moment and inverse moment priors
library(mombf)
g <- 1
thseq <- seq(-3,3,length=1000)
plot(thseq,dmom(thseq,g=g),type='l',ylab='Prior density')
lines(thseq,dimom(thseq,g=g),lty=2,col=2)
}
\keyword{ distribution }
\keyword{ htest }
