r.rcm <- function(n, Pi, rho, m)
{
	k <- length(Pi)
	stopifnot( 0 < rho && rho < 1 )
	if (length(m) == 1) m <- rep(m, n)
	Y <- rmultinom(n, prob = Pi, size = 1)
	N <- rbinom(n, prob = rho, size = m)

	## Note that in base R, size argument is not vectorized in:
	## X <- rmultinom(n, prob = Pi, size = m-N)
	X <- sapply(m-N, rmultinom, n = 1, prob = Pi)

	N.mat <- matrix(N, k, n, byrow = TRUE)
	N.mat*Y + X
}

d.rcm <- function(x, Pi, rho, m, log = FALSE)
{
	k <- length(Pi)
	Pi <- normalize(Pi) # We get strange results if Pi isn't normalized
	P <- (1 - rho) * Pi + diag(rho, k, k)
	d.multmix(x, P, Pi, m, log)
}
