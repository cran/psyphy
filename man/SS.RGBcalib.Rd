\name{SS.RGBcalib}
\alias{SS.RGBcalib}
\alias{SS.calib}

\title{  Self-Start Functions for Fitting Luminance vs Grey Level Relation on CRT displays}
\description{
  This \code{selfStart} model evaluates the parameters for describing
  the luminance vs grey level relation of the R, G and B guns of
  a CRT-like display, fitting a single exponent, gamma, for each
  of the 3 guns.  It has an initial attribute that will evaluate
  initial estimates of the parameters, \code{Blev}, \code{Br},
  \code{Bg}, \code{Bb} and \code{gamm}.
  In the case of fitting data from a single gun or for a combination of guns, as in the sum of the three for calibrating the \emph{white}, the parameter \code{k} is used for the coefficient.
  Both functions include gradient and hessian attributes.
}
\usage{
SS.calib(Blev, k, gamm, GL)
SS.RGBcalib(Blev, Br, Bg, Bb, gamm, Rgun, Ggun, Bgun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Blev}{numeric. The black level is the luminance at the 0 grey level }
  \item{k}{numeric, coefficient of one gun for fitting single gun}
  \item{Br}{numeric, coefficient of the R gun }
  \item{Bg}{numeric, coefficient of the G gun }
  \item{Bb}{numeric, coefficient of the B gun }
  \item{gamm}{numeric, the exponent, gamma, applied to the grey level }
  \item{GL}{numeric, is the grey level for the gun tested, covariate in model matrix in one gun case}
  \item{Rgun}{numeric, is a covariate in the model matrix that indicates the grey
  level for the R gun.  See the example below.}
  \item{Ggun}{numeric, is a covariate in the model matrix that indicates the grey
  level for the G gun }
  \item{Bgun}{numeric, is a covariate in the model matrix that indicates the grey
  level for the B gun }
}
\details{
  The model
  \deqn{ Lum(GL) = Blev + \beta_i * GL^\gamma }
  where i is in \{R, G, B\},
  usually provides a reasonable description of the change
  in luminance of a display gun with grey level, \code{GL}.
  This \code{selfStart}
  function estimates \eqn{\gamma} and the other parameters using the
  \code{nls} function.  It is assumed that grey level is normalized
  to the interval [0, 1].  This results in lower correlation between
  the linear coefficients of the guns, \eqn{\beta_i} , than if the
  actual bit-level is used, e.g., [0, 255], for an 8-bit graphics
  card (see the example).  
  Also, with this normalization of the data, the coefficients,
  \eqn{\beta_i}, provide estimates of the maximum luminance for
  each gun. The need for the arguments \code{Rgun}, \code{Ggun} and
  \code{Bgun} is really a kludge in order to add gradient and
  hessian information to the model.
}
\value{
  returns a numeric vector giving the estimated luminance given the parameters passed as arguments and a gradient matrix and a hessian array as attributes.}

\references{ ~put references to the literature/web site here ~ }
\author{Kenneth Knoblauch }

\seealso{\code{\link{nls}}}
\examples{

data(RGB)

#Fitting a single gun
W.nls <- nls(Lum ~ SS.calib(Blev, k, gamm, GL), data = RGB,
				subset = (Gun == "W"))
summary(W.nls)

#curvature (parameter effect) is greater when GL is 0:255
Wc.nls <- nls(Lum ~ SS.calib(Blev, k, gamm, GL*255), data = RGB,
				subset = (Gun == "W"))
MASS::rms.curv(W.nls)
MASS::rms.curv(Wc.nls)
pairs(profile(Wc.nls), absVal = FALSE)
pairs(profile(W.nls), absVal = FALSE)


#Fitting 3 guns with independent gamma's		
RGB0.nls <- nlme::nlsList(Lum ~ SS.calib(Blev, k, gamm, GL) | Gun,
		data = subset(RGB, Gun != "W"))
summary(RGB0.nls)
plot(nlme::intervals(RGB0.nls))

# Add covariates to data.frame for R, G and B grey levels
gg <- model.matrix(~-1 + Gun/GL, RGB)[ , c(5:7)]
RGB$Rgun <- gg[, 1]
RGB$Ggun <- gg[, 2]
RGB$Bgun <- gg[, 3]
RGB.nls <- nls(Lum ~ SS.RGBcalib(Blev, Br, Bg, Bb, gamm, Rgun, Ggun, Bgun), 
 			              data = RGB, subset = (Gun != "W") )
summary(RGB.nls)
confint(RGB.nls)

}
\keyword{ models }
\keyword{ nonlinear }
