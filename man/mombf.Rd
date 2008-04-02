\name{mombf}
\alias{mombf}
\alias{mombf.lm}
\alias{imombf}
\alias{imombf.lm}
\alias{zellnerbf}
\alias{zellnerbf.lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Moment and inverse moment Bayes factors for linear models. }
\description{
  \code{mombf} computes moment Bayes factors to test whether a subset of
  regression coefficients are equal to some user-specified value.
  \code{imombf} computes inverse moment Bayes factors.
  \code{zellnerbf} computes Bayes factors based on the Zellner-Siow
  prior (used to build the moment prior).
}
\usage{
mombf(lm1, coef, g, prior.mode, theta0, logbf = FALSE)
imombf(lm1,coef,g,prior.mode,nu=1,theta0,B=10^5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lm1}{ Linear model fit, as returned by \code{lm1}. }
  \item{coef}{ Vector with indexes of coefficients to be
    tested. e.g. \code{coef==c(2,3)} 
    and \code{theta0==c(0,0)} tests \code{coef(lm1)[2]=coef(lm1)[3]=0}. }
  \item{g}{ Vector with prior parameter values. See \code{dmom} and
    \code{dimom} for details. }
  \item{prior.mode}{ If specified, \code{g} is determined by calling
    \code{g2mode}. }
  \item{theta0}{ Null value for the regression coefficients. Defaults to
    0. }
  \item{logbf}{ If \code{logbf==TRUE} the natural logarithm of the Bayes
    factor is returned.}
  \item{nu}{ Prior parameter for the inverse moment prior. See
    \code{dimom} for details. Defaults to \code{nu=1}, which Cauchy-like
    tails.}
  \item{B}{Number of Monte Carlo samples to estimate the inverse moment
    Bayes factor.}
}
\details{
  These functions actually call \code{momunknown} and
  \code{imomunknown}, but they have a simpler interface.
  See \code{dmom} and \code{dimom} for details on the moment and inverse
  moment priors.
  The Zellner-Siow g-prior is given by dmvnorm(theta,theta0,n*g*V1).
}
\value{
\code{mombf} returns the moment Bayes factor to compare the model where
\code{theta!=theta0}
with the null model where \code{theta==theta0}. Large values favor the
alternative model; small values favor the null.
\code{imombf} returns
inverse moment Bayes factors.
\code{zellnerbf} returns Bayes factors based on the Zellner-Siow g-prior.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\seealso{ \code{\link{momunknown}},
  \code{\link{imomunknown}} and \code{\link{zbfunknown}} for another interface to compute Bayes
factors. \code{\link{momknown}}, \code{\link{imomknown}} and \code{\link{zbfknown}}
to compute Bayes factors assuming that the dispersion parameter
is known, and for approximate Bayes factors for
GLMs. \code{\link{mode2g}} for prior elicitation.  }
\examples{
##compute Bayes factor for Hald's data
data(hald)
lm1 <- lm(hald[,1] ~ hald[,2] + hald[,3] + hald[,4] + hald[,5])

prior.mode <- .2
V <- summary(lm1)$cov.unscaled
gmom <- mode2g.univ(prior.mode,V[2,2],nrow(hald),prior='Mom')
gimom <- mode2g.univ(prior.mode,V[2,2],nrow(hald),prior='iMom')

mombf(lm1,coef=2,g=gmom) #moment BF
imombf(lm1,coef=2,g=gimom,B=10^5) #inverse moment BF
zellnerbf(lm1,coef=2,g=1) #BF based on Zellner's g-prior

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

