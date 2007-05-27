\name{glm.lambda}
\alias{glm.lambda}
\title{ mafc Probit Fit to Psychometric Function Profiled on Upper Asymptote }
\description{
 A wrapper for \code{glm} in which the deviance for the model with binomial family and link \code{probit.lambda} is profiled as a function of \code{lambda}, the upper asymptote of the psychometric function.  
}
\usage{
glm.lambda(formula, data, NumAlt = 2, lambda = seq(0, 0.1, len = 40), 
	plot.it = FALSE, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be fit }
  \item{data}{ an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{glm} is called. }
  \item{NumAlt}{ the number of alternatives, \code{m} in the mafc experiment from which the data arise }
  \item{lambda}{ a sequence of values to profile for the upper asymptote of the psychometric function }
  \item{plot.it}{ logical indicating whether to plot the profile of the deviances as a function of \code{lambda}}
  \item{\dots}{ further arguments passed to \code{glm}}
}
\details{
  The psychometric function fit to the data is described by
  \deqn{P(x) = 1/m + (1 -  1/m - \lambda) \Phi(x)} where \eqn{m} is the number of alternatives and the lower asymptote, \eqn{1 - \lambda} is the upper asymptote and \eqn{\Phi} is the cumulative normal function.
}
\value{
  returns an object of class \sQuote{lambda} which inherits from classes \sQuote{glm} and \sQuote{lm}.  It only differs from an object of class \sQuote{glm} in including two additional components, \code{lambda}, giving the estimated minimum of the profile by fitting a quadratic to the profile and a data frame containing the profiled deviance values for each value of \code{lambda} tested.  The degrees of freedom are reduced by 1 to take into account the estimation of \code{lambda}.
  }
\references{ Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting, sampling, and goodness of fit. Percept Psychophys., 63(8), 1293--1313.

Yssaad-Fesselier, R. and Knoblauch, K. (2006) Modeling psychometric
functions in R. \emph{ Behav Res Methods.},  \bold{38(1)}, 28--41. (for examples
with \code{gnlr}).	
	}
\author{ Ken Knoblauch}
\note{ If the minimum occurs outside the interval examined, an error might occur.  In any case, re-running the function with a new range of \code{lambda} that includes the minimum should work.  if the plotted profile indicates that the fitted quadratic does not describe well the profile at the minimum, refitting with a more restricted range of \code{lambda} is recommended.
}
\seealso{ \code{\link{mafc}}, \code{\link{glm}}, \code{\link{probit.lambda}}, \code{\link{family}}}
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
set.seed(12161952)
ny <- rbinom(length(cnt), num.tr, wb(p))
nn <- num.tr - ny
phat <- ny/(ny + nn)
resp.mat <- matrix(c(ny, nn), ncol = 2)

## First with upper asymptote at 1
dd.glm <- glm(resp.mat ~ cnt, family = binomial(mafc.probit(3)))
summary(dd.glm)
dd.lam <- glm.lambda(resp.mat ~ cnt, NumAlt = 3, lambda = seq(0, 0.03, 
	len = 100), plot.it = TRUE)
summary(dd.lam)
## can refine interval, but doesn't change result much
dd.lam2 <- glm.lambda(resp.mat ~ cnt, NumAlt = 3, 
	lambda = seq(dd.lam$lambda/sqrt(2), dd.lam$lambda*sqrt(2), 
	len = 100), plot.it = TRUE)
summary(dd.lam2)
## Compare fits w/ and w/out lambda
anova(dd.glm, dd.lam2, test = "Chisq")

plot(cnt, phat, log = "x", cex = 1.5, ylim = c(0, 1))
pcnt <- seq(0.01, 0.1, len = 100)
lines(pcnt, predict(dd.glm, data.frame(cnt = pcnt),
                    type = "response"), lwd = 2)
lines(pcnt, predict(dd.lam, data.frame(cnt = pcnt),
                    type = "response"), lwd = 2, lty = 2)

}
\keyword{models}
