\name{mode2g}
\alias{mode2g}
\alias{mode2g.univ}
\alias{g2mode}
\alias{g2mode.univ}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Moment and inverse moment prior elicitation. }
\description{
  \code{mode2g} and \code{mode2g.univ} find the \code{g} value corresponding to a
  given prior mode.
  \code{g2mode} and \code{g2mode.univ} find the prior mode corresponding
  to a given \code{g} value.
}
\usage{
mode2g(prior.mode, prior = 'iMom', nu = 1, p)
mode2g.univ(prior.mode, V1, n, prior = 'iMom')
g2mode(g, n, prior = 'iMom', nu = 1, p)
g2mode.univ(g, V1, n, prior = 'iMom')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prior.mode}{ For \code{mode2g.univ}, this is the prior mode for
    |theta-theta0|/sigma, where sigma is the dispersion parameter. For
    \code{mode2g} it's the prior mode for the quadratic form
    \code{(theta-theta0)' * solve(Sigma) * (theta-theta0)/(n*g)}. }
  \item{prior}{ \code{prior=='Mom'} does computations for the moment
    prior. \code{prior=='iMom'} does computations for the inverse moment
    prior. }
  \item{nu}{ Prior parameter (only used if \code{prior=='iMom'}. }
  \item{p}{ Dimensionality of the parameter, i.e. \code{p==1} for
    univariate, \code{p==2} for bivariate and so on. }
  \item{n}{ Prior parameter. See \code{dmom} for details. }
  \item{g}{ Prior parameter. See \code{dimom} for details. }
  \item{V1}{ Prior parameter. See \code{dmom} and \code{dimom} for
    details. }
}
\details{
  See \code{dmom} and \code{dimom} for details on the meaning of the
  prior parameters.
}
\value{
  \code{mode2g} and \code{mode2g} return the value of the prior
  parameter \code{g} matching the given mode.
  \code{g2mode} and \code{g2mode.univ} return the prior mode for a given
  prior parameter \code{g}.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\seealso{ \code{\link{dmom}}, \code{\link{dimom}}, \code{\link{mombf}}, \code{\link{imombf}} }
\examples{
#find g value giving a prior mode for abs(theta/sigma) at 0.2
data(hald)
lm1 <- lm(hald[,1] ~ hald[,2] + hald[,3] + hald[,4] + hald[,5])
prior.mode <- .2
V <- summary(lm1)$cov.unscaled
gmom <- mode2g.univ(prior.mode,V[2,2],nrow(hald),prior='Mom')
gimom <- mode2g.univ(prior.mode,V[2,2],nrow(hald),prior='iMom')
gmom
gimom
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ htest }
