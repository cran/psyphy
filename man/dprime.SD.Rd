\name{dprime.SD}
\alias{dprime.SD}
\title{d' for Same-different Paradigm}
\description{
  Calulate \eqn{d'} for same-different paradigm either 
  assuming a differencing strategy or independent observations
}
\usage{
dprime.SD(H, FA, zdiff, Pcmax, method = "diff")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H}{numeric in [0, 1] corresponding to Hit rate}
  \item{FA}{numeric in [0, 1] corresponding to False alarm rate}
  \item{zdiff}{numeric. Difference of z-scores for Hit and False Alarm rates
   ( only valid for method "IO")}
  \item{Pcmax}{numeric in [0, 1].  Proportion correct for an unbiased observer,
    \code{pnorm(zdiff)} (only valid for method "IO").
  }
  \item{method}{character. Specifies the model to describe
    the observer's criterion for dividing up the decision space,
    must be either "diff" for a differencing strategy (the default)
    or "IO" for independent observations. }
}
\details{Two different strategies have been described for how the
  observer partitions the decision space in the same-different
  paradigm.  With Independent Observations, \eqn{d'} can be calculated
  either from the \code{H} and \code{FA} rates, from the difference of
  z-scores or from the probability correct of an unbiased observer.
  Only one of these three choices should be specified in the arguments.
  For the differencing strategy, only the first of these choices is valid.
}
\value{
Returns the value of \eqn{d'}
}
\references{MacMillan, N. A. and Creeman, C. D. (1991) \emph{Detection
    Theory: A User's Guide} Cambridge University Press

Green, D. M. and Swets, J. A. (1966) \emph{Signal Detection
    Theory and Psychophysics} Robert E. Krieger Publishing Company   
}
\author{Kenneth Knoblauch }

\seealso{\code{\link{dprime.mAFC}}, \code{\link{dprime.ABX}},
	\code{\link{dprime.oddity}} }
\examples{
dprime.SD(H = 0.642, F = 0.3)
dprime.SD(H = 0.75, F = 0.3, method = "IO")
dprime.SD(zdiff = qnorm(0.75) - qnorm(0.3), method = "IO")
dprime.SD(Pcmax = pnorm( (qnorm(0.75) - qnorm(0.3))/2 ), method = "IO")
}
\keyword{univar}
