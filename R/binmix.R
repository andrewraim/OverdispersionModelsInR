## Note: These don't support p or Pi varying with the observation

r.binom.mix <- function(n, p, Pi, m)
{
	J <- length(Pi)
	z <- sample(1:J, prob = Pi, replace = TRUE, size = n)

	if (class(p) == "matrix")
	{
		idx <- cbind(1:n, z)
		prob <- p[idx]
	}
	else
	{
		prob <- p[z]
	}

	rbinom(n, prob = prob, size = m)
}

d.binom.mix <- function(x, p, Pi, m, log = FALSE)
{
	J <- length(Pi)
	
	fc <- matrix(NA, length(x), J)
	for (j in 1:J)
	{
		fc[,j] <- dbinom(x, prob = p[j], size = m)
	}
	
	if (log) log(fc %*% Pi)
	else fc %*% Pi
}
