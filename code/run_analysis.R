library(tidyverse) 
source('code/utils.R')

R0 <- 2.5
Nint <- 100 
cint <- 10
cext <- 10 
pext <- 0.01
mu <- 0.01

Amean <- 2.5
Anexp <- 5

k <- cint/(cint+cext) 


inftimes <- c(-2, -2.5, -3, -4, -6, -12, -20, -25, rep(NA, Nint-8))
snpdists <- c(0, 1, 1, 2, 3, 3, 2, 5, rep(NA, Nint-8))

xi <- c(20,30,50,25,12,4,2,1)

Avals <- unlist(lapply(-inftimes, A, R0=R0, mean=2.5, nexp=5))

lambdaint <- k/Nint*sum(Avals, na.rm=TRUE)

lambdaext <- (1-k)*pext*qA(0.99, mean=Amean, nexp=Anexp)/R0

pint <- lambdaint / (lambdaint + lambdaext)

psigint <- sum(Avals * unlist(lapply(snpdists, dpois, lambda=mu)),na.rm=TRUE) / sum(Avals, na.rm=TRUE)

getpsigext <- function(xi, mu, tol=1e-6){

	xirng <- qpois(1-tol, mu) + 2

	out <- sum((xi[1:xirng] / sum(xi))*
			unlist(lapply(0:(xirng-1), dpois, lambda=mu)))

	return(out)

}

psigext <- getpsigext(xi,mu)

pintsig <- psigint*pint / (psigint*pint + psigext*(1-pint))




















xrng <- seq(from=0,to=10,by=0.01)
tibble(x=xrng, y=unlist(lapply(xrng, A, mean=Amean, nexp=Anexp))) %>% 
	ggplot(aes(x=x, y=y)) + 
		geom_line() + 
		theme_classic() 