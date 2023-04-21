A <- function(x, R0=2.5, mean=2.5, nexp=5){
	out <- R0*dgamma(x, shape=nexp, rate=nexp/mean)
	return(out)
}

qA <- function(x, mean=2.5, nexp=5){
	out <- qgamma(x, shape=nexp, rate=nexp/mean)
	return(out)
}