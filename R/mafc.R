`mafc.logit` <- 
function( m = 2 )
{
	m <- as.integer(m)
	if (m < 2) 
		stop("m must be an integer > 1")
	linkfun <- function(mu) {
		mu <- pmax(mu, 1/m +  .Machine$double.eps)
		qlogis((m * mu - 1)/(m - 1) ) }  
	linkinv <- function(eta) {
		1/m  + (m - 1)/m * .Call("logit_linkinv", eta,
            PACKAGE = "stats")
		}
	mu.eta <- function(eta) ((m -1) / m) * 
		.Call("logit_mu_eta", eta, PACKAGE = "stats")
	valideta <- function(eta) TRUE
	link <- paste("mafc.logit(", m, ")", sep = "")
	structure(list(linkfun = linkfun,
				  linkinv = linkinv,
				  mu.eta = mu.eta,
				  valideta = valideta, name = link),
				  class = "link-glm")
}

`mafc.probit` <- 
function( m = 2 )
{
	m <- as.integer(m)
	if (m < 2)
		stop("m must be an integer > 1")
	linkfun <- function(mu) {
		mu <- pmax(mu, 1/m +  .Machine$double.eps)
		qnorm((m * mu - 1)/(m - 1) ) } 
	linkinv <- function(eta) {
		1/m  + (m - 1)/m * pnorm(eta)
		}
	mu.eta <- function(eta) ((m -1) / m) * dnorm(eta)
	valideta <- function(eta) TRUE
	link <- paste("mafc.probit(", m, ")", sep = "")
	structure(list(linkfun = linkfun,
				  linkinv = linkinv,
				  mu.eta = mu.eta,
				  valideta = valideta, name = link),
				  class = "link-glm")
}
