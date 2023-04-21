# =============================================================================
# Import and define 
# =============================================================================

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

inftimes <- c(-2, -2.5, -3, -4, -6, -12, -20, -25, rep(NA, Nint-8))
snpdists <- c(0, 1, 1, 2, 3, 3, 2, 5, rep(NA, Nint-8))

xi <- c(20,30,50,25,12,4,2,1)

# =============================================================================
# Calculate probability of internal transmission
# =============================================================================

k <- cint/(cint+cext) 

Avals <- unlist(lapply(-inftimes, A, R0=R0, mean=Amean, nexp=Anexp))

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

# =============================================================================
# Plot
# =============================================================================

xrng <- seq(from=0, to=10, by=0.01) 
sdlabs <- tibble(x=inftimes+2, y=1.02, label=snpdists) %>% 
	filter(!is.na(x)) %>% 
	mutate(label=as.character(label))
fig_Aparade <- as.list(inftimes[!is.na(inftimes)]) %>% 
	imap(~ tibble(id=.y, x=xrng+.x, y=A(xrng,R0,mean=Amean,nexp=Anexp))) %>% 
	bind_rows() %>% 
	ggplot() + 
		geom_line(aes(x=x,y=y,group=factor(id)), size=0.8, alpha=0.4) + 
		geom_vline(aes(xintercept=0)) + 
		geom_hline(aes(yintercept=0)) + 
		geom_text(data=sdlabs, aes(x=x,y=y,label=label), size=2.5) + 
		theme_classic() + 
		labs(x="Time since onset of the focal infection", y="Infectiousness") +
		theme(text=element_text(size=9))
ggsave(fig_Aparade, file="figures/Aparade.pdf",width=3.5,height=3.5/1.6,units="in")
