\name{RGB}
\alias{RGB}
\docType{data}
\title{ Luminance Calibration Data from Video Projector}
\description{
  The data was obtained from the measurements of the luminance
  of the \code{R}, \code{G} and \code{B} channels individually,
  as well as the three together, W, for each of 21 grey levels,
  \code{GL} from a screen on which a video projector was displaying
  an image of a uniform field. Grey level has been normalized to
  the interval [0, 1], though originally it is specified as integers
  in [0, 255].  The measurements were obtained with a Photo Research
  650 spectro-radiometer.
}
\usage{data(RGB)}
\format{
  A data frame with 84 observations on the following 3 variables.
\describe{
    \item{\code{Lum}}{numeric vector of the measured luminance
    in candelas/meter\eqn{^2} }
    \item{\code{GL}}{The grey level normalized to the interval [0, 1]}
    \item{\code{Gun}}{factor with levels \code{R} \code{G} \code{B} \code{W}}
  }
}
\details{
}
\source{
}
\references{
}
\examples{
data(RGB)

}
\keyword{datasets}
