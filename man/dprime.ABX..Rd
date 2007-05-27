\name{dprime.ABX}
\alias{dprime.ABX}
\title{ d' for ABX Paradigm }
\description{
Calulate \eqn{d'} for ABX paradigm either 
  assuming a differencing strategy or independent observations
}
\usage{
dprime.ABX(Hits, FA, zdiff, Pc.unb, method = "diff")
}
\arguments{
  \item{Hits}{numeric in [0, 1] corresponding to Hit rate}
  \item{FA}{numeric in [0, 1] corresponding to False alarm rate}
  \item{zdiff}{numeric. Difference of z-scores for Hit and False Alarm rates
 }
  \item{Pc.unb}{numeric in [0, 1].  Proportion correct for an unbiased observer,
    \code{pnorm(zdiff)}}
  \item{method}{character. Specifies the model to describe
    the observer's criterion for dividing up the decision space,
    must be either "diff" for a differencing strategy (the default)
    or "IO" for independent observations.}
}
\details{
Two different strategies have been described for how the
  observer partitions the decision space in the ABX
  paradigm, either based on Independent Observations of
  each stimulus or on a differencing strategy.  The differecing
  strategy is the default. \eqn{d'} can be calculated
  either from the \code{H} and \code{FA} rates, from the difference of
  z-scores or from the probability correct of an unbiased observer.
}
\value{
 Returns the value of \eqn{d'}

}
\references{
  MacMillan, N. A. and Creeman, C. D. (1991) \emph{Detection
    Theory: A User's Guide} Cambridge University Press

Green, D. M. and Swets, J. A. (1966) \emph{Signal Detection
    Theory and Psychophysics} Robert E. Krieger Publishing Company  
}
\author{Kenneth Knoblauch}

\seealso{\code{\link{dprime.mAFC}}, \code{\link{dprime.SD}},
\code{\link{dprime.oddity}}}
\examples{
dprime.ABX(H = 0.75, F = 0.3)
dprime.ABX(H = 0.75, F = 0.3, method = "IO")
dprime.ABX(zdiff = qnorm(0.75) - qnorm(0.3))
dprime.ABX(Pc.unb = pnorm( (qnorm(0.75) - qnorm(0.3))/2 ))

}
\keyword{univar }
