r.dm <- function(n, Pi, rho, m)
{
	k <- length(Pi)
	Pi <- normalize(Pi)
	stopifnot( 0 < rho && rho < 1 )
	if (length(m) == 1) m <- rep(m, n)

	C <- rho^(-2) * (1 - rho^2)
	alpha <- C * Pi
	z <- rdirichlet(n, alpha)
	
	Y <- matrix(NA, k, n)
	for (i in 1:n)
	{
		Y[,i] <- rmultinom(1, prob = z[i,], size = m[i])
	}

	return(Y)
}

d.dm <- function(x, Pi, rho, m, log = FALSE)
{
	k <- length(Pi)
	Pi <- normalize(Pi)
	stopifnot( 0 < rho && rho < 1 )
	if (class(x) != "matrix") x <- as.matrix(x)
	n <- ncol(x)
	if (length(m) == 1) m <- rep(m, n)
	C <- rho^(-2) * (1 - rho^2)
	alpha <- C * Pi

	log.part <- matrix(NA, n, k)
	for (j in 1:k) {
		log.part[,j] <- lgamma(alpha[j] + x[j,]) - lgamma(alpha[j]) - lgamma(x[j,] + 1)
	}

	log.ff <- lgamma(m+1) + lgamma(C) - lgamma(C + m) + rowSums(log.part)

	if (log) return(log.ff)
	else return(exp(log.ff))
}
