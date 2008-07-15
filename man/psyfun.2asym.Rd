\name{psyfun.2asym}
\alias{psyfun.2asym}

\title{Fit Psychometric Functions and Upper and Lower Asymptotes}
\description{
 Fits psychometric functions allowing for variation of both upper and lower asymptotes.  Uses a procedure that alternates between fitting linear predictor with \code{glm} and estimating the asymptotes with \code{optim} until a minimum in -log likelihood is obtained within a tolerance. 
}
\usage{
psyfun.2asym(formula, data, link = logit.2asym, init.g = 0.01, 
	init.lam = 0.01, trace = FALSE, tol = 1e-06, 
	mxNumAlt = 50, ...)
}
\arguments{
  \item{formula}{a two sided formula specifying the response and the linear predictor }
  \item{data}{ a data frame within which the formula terms are interpreted }
  \item{link}{ a link function for the binomial family that allows specifying both upper and lower asymptotes }
  \item{init.g}{ numeric specifying the initial estimate for the lower asymptote }
  \item{init.lam}{ numeric specifying initial estimate for 1 - upper asymptote }
  \item{trace}{ logical indicating whether to show the trace of the minimization of -log likelihood }
  \item{tol}{ numeric indicating change in -log likelihood as a criterion for stopping iteration. }
  \item{mxNumAlt}{ integer indicating maximum number of alternations between \code{glm} and \code{optim} steps to perform if minimum not reached.}
  \item{\dots}{ additional arguments passed to \code{glm} }
}
\details{
  The function is a wrapper for \code{glm} for fitting psychometric functions with the equation
  \deqn{ P(x) = \gamma + (1 - \gamma - \lambda) p(x)
  	}
  where \eqn{\gamma} is the lower asymptote and \eqn{lambda} is \eqn{1 - } the upper asymptote,  and \eqn{p(x)} is the base psychometric function, varying between 0 and 1.
}
\note{The \code{cloglog.2asym} and its alias, \code{weib.2asym}, don't converge on
	occasion.  This can be observed by using the \code{trace} argument.
	One strategy is to modify the initial estimates.
	}
\value{
 list of class \sQuote{lambda} inheriting from classes \sQuote{glm} and \sQuote{lm} and containing additional components
  \item{lambda }{numeric indicating 1 - upper asymptote}
  \item{gam }{numeric indicating lower asymptote}
  \item{SElambda }{numeric indicating standard error estimate for lambda based on the Hessian of the last interation of \code{optim}.  The optimization is done on the value transformed by the function \code{plogis} and the value is stored in on this scale}
   \item{SEgam }{numeric indicating standard error estimate for gam estimated in the same fashion as \code{SElambda}}
  If a diagonal element of the Hessian is sufficiently close to 0, \code{NA} is returned.
}
\references{ Klein S. A. (2001) Measuring, estimating, and understanding the psychometric function: a commentary. \emph{Percept Psychophys.}, \bold{63(8)}, 1421--1455.

Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting, sampling, and goodness of fit. \emph{Percept Psychophys.}, \bold{63(8)}, 1293--1313. 
 }
\author{ Kenneth Knoblauch }

\seealso{ \code{\link{glm}}, \code{\link{optim}}, \code{\link{glm.lambda}}, \code{\link{mafc}} }
\examples{
#A toy example,
set.seed(12161952)
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
summary(ddprob.glm)
}


\keyword{ models }
