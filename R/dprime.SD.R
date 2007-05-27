`dprime.SD` <-
function (H, FA, zdiff, Pcmax, method = "diff") {
	if (method == "diff") {
		root2 <- sqrt(2)
		k <- root2 * qnorm(FA/2)
		est.dp <- function(dp)
			{ H - pnorm((k + dp)/root2) - pnorm((k - dp)/root2) }
		dp.res <- uniroot(est.dp, interval = c(0,5))
		dprime <- dp.res$root
	} else
	if (method == "IO")	{
		Call <- match.call() 
    if (pmatch("H", names(Call), 0) > 0) {
    	if (pmatch("FA", names(Call), 0) > 0) {
			zdiff <- qnorm(H) - qnorm(FA)
			Pcmax <- pnorm(zdiff/2)
		} else {
			zdiff <- qnorm(H) - qnorm(1-Hits)
			Pcmax <- pnorm(zdiff/2)
	} } else {
			if (pmatch("zdiff", names(Call), 0) > 0)
				{ Pcmax <- pnorm(zdiff/2)	} 
			} 
	dprime <- 2 * qnorm(0.5 * (1 + sqrt(2*Pcmax - 1)))
	} else
	{stop("method must be one of diff or IO") }
	dprime
}
