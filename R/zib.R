d.zib <- function(y, Pi, rho, m, log = FALSE)
{
	ff <- rho*(y == 0) + (1 - rho)*dbinom(y, m, Pi, log = FALSE)

	if (log) return(log(ff))
	else return(ff)
}

r.zib <- function(n, Pi, rho, m)
{
	z <- rbinom(n, size = 1, prob = rho)
	w <- rep(0,n)
	y <- rbinom(n, Pi, m)
	z*w + (1-z)*y
}
