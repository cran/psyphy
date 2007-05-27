\name{mafc}
\alias{mafc}
\alias{mafc.logit}
\alias{mafc.probit}

\title{Links for Binomial Family for m-alternative Forced-choice }
\description{ These provide links for the binomial family for fitting
  m-alternative forced-choice psychophysical functions.  }

\usage{
mafc.logit( m = 2 )
mafc.probit( m = 2 )
}

\arguments{ \item{m}{is the integer number (>1) of choices (Default to
    2AFC).  For m = 1 (Yes/No paradigm), use one of the built-in links
    for the binomial family.}  }

\details{ These
functions provide links for fitting psychometric functions arising from
an m-alternative forced-choice experiment.  The estimated coefficients
of the linear predictor influence both the location and the slope of the
psychometric function(s), but provide no means of estimating the upper
aymptote which is constrained to approach 1.  If the upper asympotote
must be estimated, it would be better to maximize directly the
likelihood, either with a function like \code{optim} or \code{gnlr} from
J. K. Lindsey's package \pkg{gnlm} (available at
\url{http://popgen.unimaas.nl/~jlindsey/rcode.html}).  }

\value{ Each
link returns a list containing functions required for relating the
response to the linear predictor in generalized linear models and the
name of the link.
\item{linkfun }{The link function}
\item{linkinv}{The inverse link function}
\item{mu.eta}{The derivative of the inverse link }
\item{valideta}{The domain over which the linear predictor is valid}
\item{link}{A name to be used for the link } }

\references{ 

Williams J, Ramaswamy D and Oulhaj A (2006) 10 Hz flicker improves
recognition memory in older people \emph{BMC Neurosci.} 2006 5;7:21
\url{http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1434755}
(for an example developed prior to this one, but for m = 2).

Klein S. A. (2001) Measuring, estimating, and understanding the
psychometric function: a commentary. \emph{Percept Psychophys.}, \bold{63(8)}, 1421--1455.

Wichmann, F. A. and Hill, N. J. (2001) The psychometric function: I.Fitting,
sampling, and goodness of fit. \emph{Percept Psychophys.},
\bold{63(8)}, 1293--1313.

Yssaad-Fesselier, R. and Knoblauch, K. (2006) Modeling psychometric
functions in R. \emph{ Behav Res Methods.},  \bold{38(1)}, 28--41. (for examples
with \code{gnlr}).

}
\author{Kenneth Knoblauch}

\seealso{
  \code{\link{family}}, \code{\link{make.link}}, \code{\link{glm}},
  \code{\link{optim}}}
}
\examples{
#A toy example,
b <- log(3.5)
g <- 1/3
d <- 0.0
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
ddlog.glm <- glm(resp.mat ~ cnt, family = binomial(mafc.logit(3)))

plot(cnt, phat, log = "x", cex = 1.5, ylim = c(0, 1))
pcnt <- seq(0.01, 0.1, len = 100)
lines(pcnt, predict(ddprob.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 2)
lines(pcnt, predict(ddlog.glm, data.frame(cnt = pcnt),
			type = "response"), lwd = 2, lty = 2)

#More interesting example with data from Yssaad-Fesselier and Knoblauch
data(ecc2)
ecc2.glm <- glm(cbind(nyes, nno) ~ Contr * Size * task, 
			family = binomial(mafc.probit(4)), data = ecc2)
summary(ecc2.glm)
ecc2$fit <- fitted(ecc2.glm)
library(lattice)
xyplot(nyes/(nyes + nno) ~ Contr | Size * task, data = ecc2,
	subscripts = TRUE, ID = with(ecc2, Size + as.numeric(task)),
	scale = list(x = list(log = TRUE), 
				 y = list(limits = c(0, 1.05))),
	xlab = "Contrast", ylab = "Proportion Correct Response",
	panel = function(x, y, subscripts, ID, ...) {
		which = unique(ID[subscripts])
		llines(x, ecc2$fit[which ==ID], col = "black", ...)
		panel.xyplot(x, y, pch = 16, ...)
		panel.abline(h = 0.25, lty = 2, ...)
		}
)

}
\keyword{ models }
