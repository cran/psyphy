\name{glm.WH}
\alias{glm.WH}
\title{ mafc Probit Fit to Psychometric Function with Upper Asymptote Less than One}
\description{
  A probit fit of a psychometric function with upper asymptote less than 1 is obtained by cycling between a fit with \code{glm} using the \code{probit.lambda} link and \code{optimize} to estimate \code{lambda}, 1 - the upper asymptotic value, until the log Likelihood changes by less than a  pre-set tolerance.
}
\usage{
glm.WH(formula, data, NumAlt = 2, lambda.init = 0.01, 
	interval = c(0, 0.05), trace = FALSE, tol = 1e-06, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ an optional data frame, list or enviroment (or object coercible by \code{\link{as.data.frame}} containing the variables in the model.  If not found in data, the variables are taken from the \code{environment(formula)}, typically the environment from \code{glm.WH} was called.}
  \item{NumAlt}{ integer indicating the number of alternatives (> 1) in the mafc-task.  (Default: 2). }
  \item{lambda.init}{ numeric, initial estimate of 1 - upper asymptote. }
  \item{interval}{ numeric vector giving interval endpoints within which to search for \code{lambda}. }
  \item{trace}{ logical, indicating whether or not to print out a trace of the iterative process. }
  \item{tol}{ numeric, tolerance for ending iterations.}
  \item{\dots}{ futher arguments passed to \code{glm}. }
}
\details{
   The psychometric function fit to the data is described by
  \deqn{P(x) = 1/m + (1 -  1/m - \lambda) \Phi(x)} where \eqn{m} is the number of alternatives and the lower asymptote, \eqn{1 - \lambda} is the upper asymptote and \eqn{\Phi} is the cumulative normal function.}
\value{
   returns an object of class \sQuote{lambda} which inherits from classes \sQuote{glm} and \sQuote{lm}.  It only differs from an object of class \sQuote{glm} in including an additional components, \code{lambda}, giving the estimated minimum of \code{lambda}.  The degrees of freedom are reduced by 1 to take into account the estimation of \code{lambda}.}
\references{Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting, sampling, and goodness of fit. Percept Psychophys., 63(8), 1293--1313.

Yssaad-Fesselier, R. and Knoblauch, K. (2006) Modeling psychometric
functions in R. \emph{ Behav Res Methods.},  \bold{38(1)}, 28--41. (for examples
with \code{gnlr}).	
 }
\author{ Ken Knoblauch }
\seealso{ \code{\link{mafc}}, \code{\link{glm}},\code{\link{glm.lambda}}, \code{\link{probit.lambda}}, \code{\link{family}} }
\examples{
b <- log(3.5)
g <- 1/4
d <- 0.04
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

tst.glm <- glm(resp.mat ~ cnt, binomial(mafc.probit(1/g)))
pcnt <- seq(0.005, 1, len = 1000)
plot(cnt, phat, log = "x", ylim = c(0, 1), xlim = c(0.005, 1),
	cex = 1.75)
lines(pcnt, predict(tst.glm, data.frame(cnt = pcnt), type = "response"), lwd = 2)
tst.lam <- glm.WH(resp.mat ~ cnt, NumAlt = 1/g, trace = TRUE)
lines(pcnt, predict(tst.lam, data.frame(cnt = pcnt), 
	type = "response"), lty = 2, lwd = 2)

}
\keyword{ models }
