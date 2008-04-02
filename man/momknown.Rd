\name{momknown}
\alias{momknown}
\alias{momunknown}
\alias{imomknown}
\alias{imomunknown}
\alias{zbfknown}
\alias{zbfunknown}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayes factors for moment, inverse moment and Zellner-Siow g-prior. }
\description{
  \code{momknown} and \code{momunknown} compute moment Bayes
  factors for linear models when \code{sigma^2} is known and unknown,
  respectively. The functions can also be used to compute approximate
  Bayes factors for generalized linear models and other settings.
  \code{imomknown}, \code{imomunknown} compute inverse
  moment Bayes factors. \code{zbfknown},
  \code{zbfunknown} compute Bayes factors based on the
  Zellner-Siow g-prior.
}
\usage{
momknown(theta1hat, V1, n, g = 1, theta0, sigma, logbf = FALSE)
momunknown(theta1hat,V1,n,nuisance.theta,g=1,theta0,ssr,logbf=FALSE)
imomknown(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,sigma,B=10^5)
imomunknown(theta1hat,V1,n,nuisance.theta,g=1,nu=1,theta0,ssr,B=10^5)
zbfknown(theta1hat,V1,n,g=1,theta0,sigma,logbf=FALSE)
zbfunknown(theta1hat,V1,n,nuisance.theta,g=1,theta0,ssr,logbf=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta1hat}{ Vector with regression coefficients estimates. }
  \item{V1}{ Matrix proportional to the covariance of
    \code{theta1hat}. For linear models, the covariance is \code{sigma^2*V1}. }
  \item{n}{ Sample size. }
  \item{nuisance.theta}{Number of nuisance regression coefficients, i.e.
    coefficients that we do not wish to test for.}
  \item{ssr}{Sum of squared residuals from a linear model call.}
  \item{g}{ Prior parameter. See \code{dmom} and \code{dimom} for details. }
  \item{theta0}{ Null value for the regression coefficients. Defaults to
    0. }
  \item{sigma}{ Dispersion parameter is \code{sigma^2}. }
  \item{logbf}{ If \code{logbf==TRUE} the natural logarithm of the Bayes
    factor is returned.}
  \item{nu}{ Prior parameter for the inverse moment prior. See
    \code{dimom} for details. Defaults to \code{nu=1}, which Cauchy-like
    tails.}
  \item{B}{Number of Monte Carlo samples to estimate the inverse moment
    Bayes factor.}  
}
\details{
  See \code{dmom} and \code{dimom} for details on the moment and inverse
  moment priors.
  The Zellner-Siow g-prior is given by dmvnorm(theta,theta0,n*g*V1).
}
\value{
\code{momknown} and \code{momunknown} return the moment Bayes factor to compare the model where
\code{theta!=theta0}
with the null model where \code{theta==theta0}. Large values favor the
alternative model; small values favor the null.
\code{imomknown} and \code{imomunknown} return
inverse moment Bayes factors.
\code{zbfknown} and \code{zbfunknown} return Bayes factors based on the Zellner-Siow g-prior.
}
\references{ See http://rosselldavid.googlepages.com for technical
  reports. }
\author{ David Rossell }
\seealso{ \code{\link{mombf}} and
  \code{\link{imombf}} for a simpler interface to compute Bayes
factors in linear regression. \code{\link{mode2g}} for prior elicitation. }
\examples{
#simulate data from probit regression
set.seed(4*2*2008)
n <- 50; theta <- c(log(2),0)
x <- matrix(NA,nrow=n,ncol=2)
x[,1] <- rnorm(n,0,1); x[,2] <- rnorm(n,.5*x[,1],1)
p <- pnorm(x[,1]*theta[1]+x[,2]+theta[2])
y <- rbinom(n,1,p)

#fit model
glm1 <- glm(y~x[,1]+x[,2],family=binomial(link = "probit"))
thetahat <- coef(glm1)
V <- summary(glm1)$cov.scaled

#compute Bayes factors to test whether x[,1] can be dropped from the model
g <- .5
bfmom.1 <- momknown(thetahat[2],V[2,2],n=n,g=g,sigma=1)
bfimom.1 <- imomknown(thetahat[2],V[2,2],n=n,nuisance.theta=2,g=g,sigma=1)
bfmom.1
bfimom.1

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ htest }

