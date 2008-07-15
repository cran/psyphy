\name{logit.2asym}
\alias{logit.2asym}
\alias{probit.2asym}
\alias{cauchit.2asym}
\alias{cloglog.2asym}
\alias{weib.2asym}
\title{ Links for Binomial Family with Variable Upper/Lower Asymptotes}
\description{
 These functions provide links for the binamial family so that psychometric functions can be fit with \emph{both} the upper and lower asymptotes different from 1 and 0, respectively.  
}
\usage{
logit.2asym(g, lam)
probit.2asym(g, lam)
cauchit.2asym(g, lam)
cloglog.2asym(g, lam)
weib.2asym( ... )
}
\arguments{
  \item{g}{ numeric in the range (0, 1), normally <= 0.5, however, which specifies the lower asymptote of the psychometric function. }
  \item{lam}{ numeric in the range (0, 1), specifying 1 - the upper asymptote of the psychometric function. }
  \item{...}{ used just to pass along the formals of \code{cloglog.2asym} as arguments to \code{weib.2asym}. }
}
\details{
  These links are used to specify psychometric functions with the form
  \deqn{ P(x) = \gamma + (1 - \gamma - \lambda) p(x)
  	}
  where \eqn{\gamma} is the lower asymptote and \eqn{lambda} is \eqn{1 - } the upper asymptote,  and \eqn{p(x)} is the base psychometric function, varying between 0 and 1.
}
\value{
  Each
link returns a list containing functions required for relating the
response to the linear predictor in generalized linear models and the
name of the link.
\item{linkfun }{The link function}
\item{linkinv}{The inverse link function}
\item{mu.eta}{The derivative of the inverse link }
\item{valideta}{The domain over which the linear predictor is valid}
\item{link}{A name to be used for the link }}
\references{ 
Klein S. A. (2001) Measuring, estimating, and understanding the
psychometric function: a commentary. \emph{Percept Psychophys.}, \bold{63(8)}, 1421--1455.

Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting, sampling, and goodness of fit. \emph{Percept Psychophys.}, \bold{63(8)}, 1293--1313. }
\author{Kenneth Knoblauch}
\seealso{\code{\link{glm}}, \code{\link{glm}} \code{\link{make.link}}, \code{\link{psyfun.2asym}}}
\examples{
#A toy example,
b <- 3
g <- 0.05 # simulated false alarm rate
d <- 0.03
a <- 0.04
p <- c(a, b, g, d)
num.tr <- 160
cnt <- 10^seq(-2, -1, length = 6) # contrast levels

#simulated Weibull-Quick observer responses
truep <- g + (1 - g - d) * pweibull(cnt, b, a)
ny <- rbinom(length(cnt), num.tr, truep)
nn <- num.tr - ny
phat <- ny/(ny + nn)
resp.mat <- matrix(c(ny, nn), ncol = 2)

ddprob.glm <- psyfun.2asym(resp.mat ~ cnt, link = probit.2asym)
ddlog.glm <- psyfun.2asym(resp.mat ~ cnt, link = logit.2asym)
# Can fit a Weibull function, but use log contrast as variable
ddweib.glm <- psyfun.2asym(resp.mat ~ log(cnt), link = weib.2asym) 
ddcau.glm <- psyfun.2asym(resp.mat ~ cnt, link = cauchit.2asym)

plot(cnt, phat, log = "x", cex = 1.5, ylim = c(0, 1))
pcnt <- seq(0.01, 0.1, len = 100)
lines(pcnt, predict(ddprob.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 5)
lines(pcnt, predict(ddlog.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 2, lty = 2, col = "blue")
lines(pcnt, predict(ddweib.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 3, col = "grey")
lines(pcnt, predict(ddcau.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 3, col = "grey", lty = 2)

}
\keyword{ models }
