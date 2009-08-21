\name{ecc2}
\alias{ecc2}
\docType{data}
\title{4-afc Detection and Identification of Letters}
\description{
  Letter detection and identification at 2 degrees eccentricity in
  the visual field.  On each trial, one of four letters (b, d, p, q) were presented in one of four positions (superior, inferior, left, right) in the visual field.  In a given session, the letter
  height was fixed.  Six contrast levels were tested in each
  session.  The data indicate the proportion of correctly
  identified positions, referred to here as detection, and
  the proportion of correctly identified letters, conditional
  on correct identification.
}
\usage{data(ecc2)}
\format{
  A data frame with 48 observations on the following 5 variables.
  \describe{
    \item{\code{Contr}}{numeric. The contrast of the stimulus, defined as
    	Weberian contrast.}
    \item{\code{task}}{a factor with levels \code{DET} \code{ID}
    	indicating the two tasks, detection and identification.}
    \item{\code{Size}}{a numeric vector indicating the letter height}
    \item{\code{Correct}}{an integer vector indicating the number of correct responses (\code{DET} or \code{ID}).}
    \item{\code{Incorrect}}{an integer vector, indicating the number
    	of incorrect responses.}
  }
}

\references{
   Yssaad-Fesselier, R. and Knoblauch, K. (2006) Modeling psychometric
functions in R. \emph{ Behav Res Methods.},  \bold{38(1)}, 28--41.
}
\examples{
data(ecc2)
library(lattice)
xyplot(Correct/(Correct + Incorrect) ~ Contr | Size * task, ecc2, 
		type = "b", scale = list(x = list(log = TRUE), 
		y = list(limits = c(0, 1.05))),
		xlab = "Contrast", ylab = "Proportion Correct Response",
		panel = function(x, y, ...) {
			panel.xyplot(x, y, ...)
			panel.abline(h = 0.25, lty = 2)
			})

}
\keyword{datasets}
