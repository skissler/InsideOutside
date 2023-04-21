library(tidyverse) 
source('code/utils.R')

R0 <- 2.5
Nint <- 100 
cint <- 10
cext <- 10 
pext <- 0.01

Amean <- 2.5
Anexp <- 5

k <- cint/(cint+cext) 


inftimes <- c(-2, -2.5, -3, -4, -6, -12, -20, -25, rep(NA, Nint-8))


lambdaint <- k/Nint*sum(unlist(lapply(
	-inftimes, A, R0=R0, mean=2.5, nexp=5)),
	na.rm=TRUE)

lambdaext <- (1-k)*pext*qA(0.99, mean=Amean, nexp=Anexp)/R0


pint <- lambdaint / (lambdaint + lambdaext)


xrng <- seq(from=0,to=10,by=0.01)
tibble(x=xrng, y=unlist(lapply(xrng, A, mean=Amean, nexp=Anexp))) %>% 
	ggplot(aes(x=x, y=y)) + 
		geom_line() + 
		theme_classic() 