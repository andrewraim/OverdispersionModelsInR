# --- My calculation for the expected value ---
# (1 - Pi) * mu / (1 - 1/(1 + kappa * mu)^(1/kappa))

# mu <- 5
# kappa <- 3
# Pi <- 0.05

# --- Using the density --
# y <- 0:500
# ff <- d.hurdle.nb(y, mu, kappa, Pi)
# E <- sum(y * ff)

# --- Using the drawing function ---
# y <- r.hurdle.nb(1000000, mu, kappa, Pi)
# mean(y)

r.hurdle.nb <- function(n, mu, kappa, Pi)
{
	if (length(kappa) == 1) kappa <- rep(kappa, n)
	if (length(mu) == 1) mu <- rep(mu, n)

    z <- rbinom(n, 1, Pi)

	# To draw from truncated NB, condition on drawing x > 0
	x <- rnbinom(n, mu = mu, size = 1/kappa)
	idx <- which(x == 0)
	while (length(idx) > 0) {
		x[idx] <- rnbinom(length(idx), mu = mu[idx], size = 1/kappa[idx])
		idx <- which(x == 0)
	}

	(z==1)*0 + (z==0)*x
}

d.hurdle.nb <- function(x, mu, kappa, Pi, log = FALSE)
{
	fz <- (x == 0)
	f1 <- dnbinom(x, mu = mu, size = 1/kappa, log = FALSE) * (x != 0)
	f0 <- dnbinom(0, mu = mu, size = 1/kappa, log = FALSE)
	ff <- Pi*fz + (1-Pi)*f1/(1-f0)
	if (log) return(log(ff))
	else return(ff)
}
