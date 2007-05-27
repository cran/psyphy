\name{probit.lambda}
\alias{probit.lambda}
\title{ mafc Probit Link for Binomial Family with Upper Asymptote < 1 }
\description{
This function provides a link for the binomial family for fitting m-alternative forced-choice, with a probit link and with the upper asymptote permitted to be less than 1.}
\usage{
probit.lambda(m = 2, lambda = 0)
}
\arguments{
  \item{m}{ is the integer number (>1) of choices (Default to 2AFC).  }
  \item{lambda}{ number in [0, 1] indicating 1 minus the upper asymptotic value of the psychometric function. }
}
\details{
  This function provides a link for fitting psychometric functions arising from an m-alternative forced-choice experiment using a probit link and allowing that the upper aymptote is less than 1.  The psychometric function fit to the data is described by
  \deqn{P(x) = 1/m + (1 -  1/m - \lambda) \Phi(x)} where \eqn{m} is the number of alternatives and the lower asymptote, \eqn{1 - \lambda} is the upper asymptote and \eqn{\Phi} is the cumulative normal function.
}
\value{
The link returns a list containing functions required for relating the response to the linear predictor in generalized linear models and the name of the link.

  \item{linkfun }{The link function}
  \item{linkinv }{DTHe inverse link function}
  \item{mu.eta }{The derivative of the inverse link function}
  \item{valideta }{The domain over which the linear predictor is valid}
  \item{link }{A name to be used for the link}
}
\references{Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting, sampling, and goodness of fit. Percept Psychophys., 63(8), 1293--1313.}
\author{Ken Knoblauch}
\note{Due to the difficulty of the task, subject error or incorrectly recorded data, psychophysical data may reveal less than perfect performance when stimulus differences are readily visible.  When this occurs, letting the upper asymptote be less than 1 often results in a better fit to the data and a less-biased estimate of the steepness of the curve (see example below). 
}
\seealso{\code{\link{mafc}}, \code{\link{glm}}, \code{\link{glm.lambda}}, \code{\link{family}}, \code{\link{make.link}} }
\examples{
b <- log(3.5)
g <- 1/3
d <- 0.025
a <- 0.04
p <- c(a, b, g, d)
num.tr <- 160
cnt <- 10^seq(-2, -1, length = 6) # contrast levels

#simulated observer responses
wb <- function(p) {
                   p[3] + (1 - p[3] - p[4]) *
                        (1 - exp(-((cnt/p[1])^exp(p[2]))))
                   }
ny <- rbinom(length(cnt), num.tr, wb(p))
nn <- num.tr - ny
phat <- ny/(ny + nn)
resp.mat <- matrix(c(ny, nn), ncol = 2)

ddprob.glm <- glm(resp.mat ~ cnt, family = binomial(mafc.probit(3)))
ddprob.lam <- glm(resp.mat ~ cnt, family = binomial(probit.lambda(3, 0.025)))
AIC(ddprob.glm, ddprob.lam)

plot(cnt, phat, log = "x", cex = 1.5, ylim = c(0, 1))
pcnt <- seq(0.01, 0.1, len = 100)
lines(pcnt, predict(ddprob.glm, data.frame(cnt = pcnt),
                        type = "response"), lwd = 2)
lines(pcnt, predict(ddprob.lam, data.frame(cnt = pcnt),
                        type = "response"), lwd = 2, lty = 2)
}
\keyword{models}
