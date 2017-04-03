# --- My calculation for the expected value ---
# (1 - Pi) * mu / (1 - exp(-mu))

# mu <- 5
# Pi <- 0.05

# --- Using the density --
# y <- 0:500
# ff <- d.hurdle.pois(y, mu, Pi)
# E <- sum(y * ff)

# --- Using the drawing function ---
# y <- r.hurdle.pois(1000000, mu, Pi)
# mean(y)

r.hurdle.pois <- function(n, mu, Pi)
{
	if (length(mu) == 1) mu <- rep(mu, n)

    z <- rbinom(n, 1, Pi)

	# To draw from truncated NB, condition on drawing x > 0
	x <- rpois(n, mu)
	idx <- which(x == 0)
	while (length(idx) > 0) {
		x[idx] <- rpois(length(idx), mu[idx])
		idx <- which(x == 0)
	}

	(z==1)*0 + (z==0)*x
}

d.hurdle.pois <- function(x, mu, Pi, log = FALSE)
{
	fz <- (x == 0)
	f1 <- dpois(x, mu, log = FALSE) * (x != 0)
	f0 <- dpois(0, mu, log = FALSE)
	ff <- Pi*fz + (1-Pi)*f1/(1-f0)
	if (log) return(log(ff))
	else return(ff)
}
