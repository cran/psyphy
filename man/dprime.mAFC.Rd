\name{dprime.mAFC}
\alias{dprime.mAFC}
\title{ d' for m-alternative Forced-choice }
\description{
Calculate the value of \eqn{d'} for an m-alternative forced choice paradigm
}
\usage{
dprime.mAFC(Pc, m)
}
\arguments{
  \item{Pc}{ The proportion of correct responses based on either the Hit
  rate or based on an unbiased observer}
  \item{m}{ The number of alternative choices, an integer > 1. }
}
\details{
  The probability of a correct response in m-alternative forced-choice,
  assuming independence, is based on the product of the likelihoods
  of the signal alternative generating the strongest response and
  the m - 1 noise alternatives generating responses less than this
  value (Green and Dai, 1991).
  For a Yes-No paradigm, the sensitivity is calculated
  more simply as \deqn{d' = \code{qnorm}(H) - \code{qnorm}(F)}
  where \eqn{H} and \eqn{F} are the Hit and False Alarm rates,
  respectively.
}
\value{
 Returns the value of \eqn{d'}
}
\references{
  Green, D. M. and Dai, H. (1991) Probability of being correct with
  1 of M orthogonal signals.  \emph{Perception & Psychophysics},
  \bold{49}, 100--101.
  
  Green, D. M. and Swets, J. A. (1966) \emph{Signal Detection
    Theory and Psychophysics} Robert E. Krieger Publishing Company  }
\author{ Kenneth Knoblauch }
\note{  
  Currently is only valid for \eqn{d'} in the interval [-10, 10]
  which should be well outside the range of sensory differences that
  this paradigm is used to investigate.
}
\seealso{ See Also \code{\link{dprime.ABX}}, \code{\link{dprime.SD}},
	\code{\link{dprime.oddity}}}
\examples{

dprime.mAFC(0.8, 4)
}
\keyword{univar}

