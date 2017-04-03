r.zip <- function(n, mu, Pi)
{
	z <- rbinom(n, 1, Pi)
	x <- rpois(n, mu)
	(z==1)*0 + (z==0)*x
}

d.zip <- function(x, mu, Pi, log = FALSE)
{
    fz <- (x == 0)
    f1 <- dpois(x, lambda = mu, log = FALSE)
    ff <- Pi*fz + (1-Pi)*f1
    if (log) return(log(ff))
    else return(ff)
}
