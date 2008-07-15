psyfun.2asym <- function(formula, data, link = logit.2asym,
	init.g = 0.01, init.lam = 0.01,
	trace = FALSE, tol = 1e-6, mxNumAlt = 50, ...) {
	p.l <- function(p) {
        link(p[1], p[2])
        }
	if (missing(data)) 
        data <- environment(formula)
    est.glm <- glm(formula, 
    	family = binomial(link(g = init.g, lam = init.lam)), 
    	data = data, ...)
    ll <- function(p, x) {
    	p <- plogis(p)
        rr <- x$model[[1]]
        -sum(rr[, 1] * log(p.l(p)$linkinv(x$linear.predictors)) + 
            rr[, 2] * log(1 - p.l(p)$linkinv(x$linear.predictors)))
    }
    dlogcur <- dd <- -as.vector(logLik(est.glm))
    new.glm <- est.glm
    n <- 0
    p <- c(init.g, init.lam)
	while (dd > tol) {
        n <- n + 1
        p <- qlogis(p)
        p <- optim(p, ll, x = new.glm)$par
        p <- plogis(p)
        new.glm <- glm(formula, 
        	family = binomial(link(g = p[1], lam = p[2])), 
        	data = data, ...)
        dd <- abs(-as.vector(logLik(new.glm)) - dlogcur)/dlogcur
        dlogcur <- -as.vector(logLik(new.glm))
        if (trace) 
            print(data.frame(n = n, logLik = dd, lambda = p[2], gamma = p[1]))
        if (n > mxNumAlt) {
        	print("Number of iterations exceeded without finding best fit. \n")
        	break}    
    }
	new.glm$lambda <- p[2]
	new.glm$gam <- p[1]
	cat("lambda = \t", p[2], "\t", "gamma = ", p[1], "\n")
	new.glm$df.residual <- new.glm$df.residual - 2
	new.glm$call[[3]][[2]][[1]] <- as.name(substitute(link))
	class(new.glm) <- c("lambda", "glm", "lm")
	new.glm
}