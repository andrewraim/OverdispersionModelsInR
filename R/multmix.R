r.multmix <- function(n, P, Pi, m)
{
	v <- length(Pi)
	k <- nrow(P)
	z <- sample(1:v, prob = Pi, replace = TRUE, size = n)
	
	if (length(m) == 1) m <- rep(m,n)
	y <- matrix(NA, nrow = k, ncol = n)

	for (i in 1:n) {
		y[,i] <- rmultinom(1, prob = P[,z[i]], size = m[i])
	}

	return(y)
}

## A vectorized version of the Multinomial density. Is effective when there
## are many observations and not too many categories
dmultinom <- function(x, size, prob, log = FALSE)
{
	if (class(x) != "matrix") x <- as.matrix(x)
	k <- nrow(x)
	n <- ncol(x)

	if (length(size) == 1) m <- rep(size,n)
	else m <- size

	fj <- matrix(NA, n, k)
	for (j in 1:k) {
		fj[,j] <- -lgamma(x[j,]+1) + x[j,] * log(prob[j])
	}

	log.ff <- lgamma(m+1) + rowSums(fj)

	if (log) return(log.ff)
	else return(exp(log.ff))
}

d.multmix <- function(x, P, Pi, m, log = FALSE)
{
	Pi <- normalize(Pi)
	if (class(x) != "matrix") x <- as.matrix(x)
	J <- length(Pi)
	n <- ncol(x)
	if (length(m) == 1) m <- rep(m,n)
	
	ff <- matrix(NA, n, J)
	for(j in 1:J) {
		ff[,j] <- dmultinom(x, m, P[,j])
		#for (i in 1:n)
		#{
		#	ff[i,j] <- dmultinom(x[,i], m[i], P[,j])
		#}
	}
	
	res <- as.numeric(ff %*% Pi)

	if (log) return(log(res))
	else return(res)
}
