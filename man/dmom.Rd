\name{dmom}
\alias{dmom}
\alias{dimom}
\alias{pmom}
\alias{pimom}
\alias{qmom}
\alias{qimom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Moment prior and inverse moment prior. }
\description{
  \code{dmom} and \code{dimom} return the density for the moment and inverse moment
  priors, in the univariate and multivariate setting.
  \code{pmom} and \code{pimom} return the distribution function for the univariate
  moment and inverse moment priors.
  \code{qmom} and \code{qimom} return the quantiles for the univariate
  moment and inverse moment priors.
}
\usage{
dmom(x, V1 = 1, g = 1, n = 1, theta0, baseDensity = "normal", nu = 3)
dimom(x, V1 = 1, g = 1, n = 1, nu = 1, theta0, logscale = FALSE)
pmom(q, V1 = 1, g = 1, n = 1)
pimom(q, V1 = 1, g = 1, n = 1, nu = 1)
qmom(p, V1 = 1, g = 1, n = 1)
qimom(p, V1 = 1, g = 1, n = 1, nu = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ In the univariate setting, \code{x} is a vector with the
    values at which to evaluate the density. In the multivariate setting
    it is a matrix with an observation in each row.}
  \item{q}{ Vector of quantiles.}
  \item{p}{ Vector of probabilities.}
  \item{theta0}{ Location parameter. Defaults to 0. }
  \item{V1}{ Scale matrix. Defaults to 1 in univariate setting and
    the identity matrix in the multivariate setting. }
  \item{g}{ Prior parameter. See details. }
  \item{n}{ Prior parameter. See details. }
  \item{baseDensity}{For \code{baseDensity=='normal'} a normal MOM prior
  is used, for \code{baseDensity=='t'} a T MOM prior with \code{nu}
  degrees of freedom is used.}
  \item{nu}{ Prior parameter indicating the degrees of freedom for the T
    MOM and iMOM prior densities. The
    tails of the inverse moment prior are proportional to the tails of a
    multivariate T with \code{nu} degrees of freedom.}
  \item{logscale}{ For \code{logscale==TRUE}, \code{dimom} returns the
    natural log of the prior density.} 
}
\details{
 Define the quadratic form q(theta)= (theta-theta0)' *
   solve(V1) * (theta-theta0) / (n*g).
   The normal moment prior density is proportional to
   q(theta)*dmvnorm(theta,theta0,n*g*V1).
   The T moment prior is proportional to
   q(theta)*dmvt(theta,theta0,n*g*V1,df=nu).
   The inverse moment prior density is proportional to
   \code{q(theta)^(-(nu+d)/2) * exp(-1/q(theta))}.
   pmom, pimom and qimom use closed-form expressions, while qmom uses
   nlminb to find quantiles numerically.
}
\value{
  \code{dmom} returns the value of the moment prior density.
  \code{dimom} returns the value of the inverse moment prior density.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\seealso{ \code{\link{g2mode}} to find the
  prior mode corresponding to a given \code{g}. \code{\link{mode2g}}
  to find the g value corresponding to a given prior mode.}
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
