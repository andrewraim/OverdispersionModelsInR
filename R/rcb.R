r.rcb <- function(n, Pi, rho, m)
{
	stopifnot( 0 < rho && rho < 1 )
	Y <- rbinom(n, prob = Pi, size = 1)
	N <- rbinom(n, prob = rho, size = m)
	X <- rbinom(n, prob = Pi, size = m-N)
	N*Y + X
}

d.rcb <- function(x, Pi, rho, m, log = FALSE)
{
	stopifnot( 0 < rho && rho < 1 )
	fc1 <- dbinom(x, prob = (1 - rho) * Pi + rho, size = m)
	fc2 <- dbinom(x, prob = (1 - rho) * Pi, size = m)
	ff <- Pi * fc1 + (1-Pi) * fc2
	
	if (log) return(log(ff))
	else return(ff)
}
