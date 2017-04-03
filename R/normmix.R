r.normmix <- function(n, mu, sigma, Pi)
{
	J <- length(Pi)
	z <- sample(1:J, prob = Pi, replace = TRUE, size = n)
	y <- rnorm(n, mu[z], sigma[z])
	return(y)
}

d.normmix <- function(x, mu, sigma, Pi, log = FALSE)
{
	Pi <- normalize(Pi)
	J <- length(Pi)
	n <- length(x)
	
	ff <- matrix(NA, n, J)
	for(j in 1:J) {
		ff[,j] <- dnorm(x, mu[j], sigma[j])
	}

	res <- as.numeric(ff %*% Pi)

	if (log) return(log(res))
	else return(res)
}
