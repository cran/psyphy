\name{dprime.oddity}
\alias{dprime.oddity}
\title{d' for 3 Stimulus Oddity Paradigm }
\description{
  Calculate \eqn{d'} for a 3 stimulus (triangular) paradigm.  Two of the
  stimuli are the same and the observer must designate the stimulus that
  is different.
}
\usage{
dprime.oddity(Pc.tri)
}
\arguments{
  \item{Pc.tri}{numeric in (1/3, 1).  The proportion of correct
    responses for an unbiased observer.
}
}
\details{

}
\value{
 Returns the value of \eqn{d'}
}
\references{
  Frijters, G. S., Kooistra, A. and Verijken, P. F. G. (1980) Tables of \eqn{d'}
  for the triangular method and the 3-AFC signal detection procedure.
  \emph{Perception & Psychophysics}, \bold{27}, 176--178.
  
MacMillan, N. A. and Creeman, C. D. (1991) \emph{Detection
    Theory: A User's Guide} Cambridge University Press

Green, D. M. and Swets, J. A. (1966) \emph{Signal Detection
    Theory and Psychophysics} Robert E. Krieger Publishing Company  
}
\author{Kenneth Knoblauch}

\seealso{\code{\link{dprime.mAFC}}, \code{\link{dprime.SD}},
          \code{\link{dprime.ABX}}
}
\examples{
dprime.oddity(0.8)
}
\keyword{univar}
