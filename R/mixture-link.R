## Andrew Raim
## Mixture Link Binomial code
## This was adapted from my thesis code, and still needs some cleanup. There are
## currently several approaches used to evaluate the density: a method based on
## Imhof (1961), a beta approximation evaluated via quadrature or the integrate
## function, and a basic monte carlo (the least efficient)
##
## Imhof works well when J=2, but becomes extremely time consuming when J>2.
## When J>2, the beta approximation is more efficient, and also seems to give
## an accurate result.

# ------------------- Draw from polytope using Dirichlet -------------------
r.dirichlet.polytope.x <- function(n, x, Beta, Pi, alpha)
{
	p <- plogis(t(x) %*% Beta)
	r.dirichlet.polytope(n, p, Pi, alpha)
}

r.dirichlet.polytope <- function(n, p, Pi, alpha)
{
	V <- find.vertices(p, Pi)
	lambda <- rdirichlet(n, alpha)
	y <- V %*% t(lambda)
	return(t(y))
}

# ------------------- Draw from Mixture Link -------------------
r.mixture.link.x <- function(X, Beta, Pi, kappa, m)
{
	n <- nrow(X)
	if (length(m) == 1) m <- rep(m, n)
	p <- plogis(X %*% Beta)

	# Note: This does not work for varying p
	# r.mixture.link(n, p, Pi, kappa, m)

	y <- numeric(n)
	for (i in 1:n)
	{
		y[i] <- r.mixture.link(1, p[i], Pi, kappa, m[i])
	}
	
	return(y)
}

r.mixture.link <- function(n, p, Pi, kappa, m)
{
	stopifnot(length(p) == 1)

	# Warning: This function currently does not work if p is a vector
	if (length(m) == 1) m <- rep(m, n)

	J <- length(Pi)
	V <- find.vertices(p, Pi)
	k <- ncol(V)
	alpha <- rep(kappa, k)
	lambda <- rdirichlet(n, alpha)
	mu <- t(V %*% t(lambda))
	r.binom.mix(n, p = mu, Pi = Pi, m = m)
}

# ------------------- Find Vertices -------------------
find.vertices.x <- function(x, Beta, Pi, tol = 1e-8)
{
	p <- plogis(t(x) %*% Beta)
	find.vertices(p, Pi, tol)
}

find.vertices <- function(p, Pi, tol = 1e-8)
{
	vert <- .Call("find_vertices", p, Pi, tol)
	return(vert)
}

# ------------------- Dist'n of a Linear Comb. of Dirichlet -------------------
p.lincomb.dirichlet.kernel.one <- function(u, x, a, alpha)
{
	h <- 2 * alpha
	Theta <- 1/2 * sum(h * atan((a-x)*u))
	Rho <- prod( (1 + (a-x)^2 * u^2)^(h/4) )
	ff <- sin(Theta) / (u*Rho)

	return(ff)
}
p.lincomb.dirichlet.kernel <- Vectorize(p.lincomb.dirichlet.kernel.one, vectorize.args = c("u"))

p.lincomb.dirichlet.one <- function(x, a, alpha, method = "imhof")
{
	if (method == "integrate")
	{
		II <- integrate(p.lincomb.dirichlet.kernel, x = x, a = a,
			alpha = alpha, lower = 0, upper = Inf,
			abs.tol = 1e-10)
			#rel.tol = 55*.Machine$double.eps)
		ff <- 1/2 - 1/pi * II$value
	}
	else if (method == "imhof")
	{
		ss <- imhof(q = 0, lambda = a-x, h = 2*alpha)
		ff <- 1 - ss$Qq
	}
	else
	{
		stop("Unknown method specified in p.lincomb.dirichlet.one")
	}

	return(ff)
}
p.lincomb.dirichlet <- Vectorize(p.lincomb.dirichlet.one, vectorize.args = c("x"))

# A more general form than k=2 case from appendix of Provost and Cheung (2000)
d.lincomb.dirichlet.k2 <- function(x, a, alpha, log = FALSE)
{
	# Assuming below that the dist'n was derived s.t. a[2] < a[1]
	# If that's not the case, reverse here
	if (a[1] < a[2])
	{
		a <- rev(a)
		alpha <- rev(alpha)
	}

	log.ff <- rep(-Inf, length(x))
	idx <- which( a[2] <= x & x <= a[1] )

	log.ff[idx] <- -lbeta(alpha[2], alpha[1]) +
		(alpha[2] - 1) * log(a[1] - x[idx]) +
		(alpha[1] - 1) * log(x[idx] - a[2]) +
		(1 - alpha[1] - alpha[2]) * log(a[1] - a[2])

	if (log) return(log.ff)
	else return(exp(log.ff))
}

d.lincomb.dirichlet <- function(x, a, alpha, eps = 1e-6, log = FALSE, method = "imhof")
{
	k <- length(a)

	if (k == 2)
	{
		ff <- d.lincomb.dirichlet.k2(x, a, alpha)
	}
	else
	{
		# Use simple finite difference to get density from CDF
		# Note that this can lead to numerically negative densities, so
		# we apply a correction at the end
		f1 <- p.lincomb.dirichlet(x + eps, a, alpha, method)
		f2 <- p.lincomb.dirichlet(x, a, alpha, method)
		ff <- pmax((f1 - f2) / eps, 0)
	}

	if (log) return(log(ff))
	else return(ff)	
}

# ------------------------- Mixture Link Dist'n -------------------------

d.mixture.link.one.x <- function(y, x, m, Beta, Pi, kappa, log = FALSE,
	method = "imhof")
{
	p <- plogis( t(x) %*% Beta )
	d.mixture.link.one(y, m, p, Pi, kappa, log, method)
}

# Use integrate function to compute the beta approx to the density
# Seems much more efficient than simple quadrature, but sometimes fails to converge
d.mixture.link.one.betaapprox <- function(y, m, p, Pi, kappa, log)
{
	V <- find.vertices(p, Pi)
	d.mixture.link.one.betaapprox.raw(y, m, V, Pi, kappa, log)
}

d.mixture.link.one.betaapprox.kernel <- function(w, a, b, lo, hi, y, m, const = NULL)
{
	log.m.choose.y <- ifelse(is.null(const), lgamma(m+1) - lgamma(y+1) - lgamma(m-y+1), const)
	log.integrand <- log.m.choose.y + y*log(w) + (m-y)*log(1-w)
	log.density <- -log(hi - lo) + dbeta((w - lo) / (hi - lo), a, b, log = TRUE)
	exp(log.integrand + log.density)
}

d.mixture.link.one.betaapprox.raw <- function(y, m, V, Pi, kappa, log)
{
	stopifnot (length(kappa) == 1)

	J <- nrow(V)
	k <- ncol(V)

	# E() and Var() of lin comb of Dirichlet	
	xi <- rowMeans(V)
	tau.sq <- numeric(J)
	ff <- numeric(J)

	for (j in 1:J)
	{
		num <- k * t(V[j,]) %*% V[j,] - (k*xi[j])^2
		denom <- k^2 * (1 + k * kappa)
		tau.sq[j] <- num / denom
	}

	lo <- apply(V, 1, min)
	hi <- apply(V, 1, max)

	# Equate beta moments to xi and tau.sq for j = 1, ..., J
	a <- 1/tau.sq * (xi - lo)^2 * (hi-xi)/(hi-lo) - (xi-lo)/(hi-lo)
	b <- a * (hi-xi)/(xi-lo)

	log.m.choose.y <- lgamma(m+1) - lgamma(y+1) - lgamma(m-y+1)
	is.degenerate <- ( sapply(hi-lo, all.equal, 0) == TRUE )

	for (j in 1:J)
	{
		# If hi[j] == lo[j], integration is over a singleton set
		if (is.degenerate[j])
		{
			ff[j] <- exp( log.m.choose.y + y*log(lo[j]) + (m-y)*log(1-lo[j]) )
		}
		else
		{
			II <- integrate(d.mixture.link.one.betaapprox.kernel, lower = lo[j],
				upper = hi[j], a = a[j], b = b[j], lo = lo[j], hi = hi[j], y = y,
				m = m, const = log.m.choose.y)
			ff[j] <- II$value
		}
	}

	res <- sum(Pi * ff)

	if (log) return(log(res))
	else return(res)
}

d.mixture.link.one.betaapprox <- function(y, m, p, Pi, kappa, log)
{
	V <- find.vertices(p, Pi)
	d.mixture.link.one.betaapprox.raw(y, m, V, Pi, kappa, log)
}
d.mixture.link.betaapprox <- Vectorize(d.mixture.link.one.betaapprox,
	vectorize.args = c("y", "m"))

d.mixture.link.one <- function(y, m, p, Pi, kappa, log = FALSE,
	method = "imhof", mc.reps = 10000)
{
	J <- length(Pi)
	V <- find.vertices(p, Pi)
	k <- ncol(V)
	stopifnot (length(kappa) == 1)
	alpha <- rep(kappa, k)

	if (k == 1)
	{
		# This is a degenerate situation. There's only one point in the set A,
		# so the integral goes away and we just evaluate at this point
		mu <- V
		fc <- exp(y * log(mu) + (m-y) * log(1-mu))
		ff <- choose(m,y) * sum(Pi * fc)
	}
	else if (method == "betaapprox")
	{
		ff <- d.mixture.link.one.betaapprox.raw(y, m, V, Pi, kappa, log = FALSE)
	}
	else if (method %in% c("integrate", "imhof"))
	{
		g <- function(x, a, alpha)
		{
			# ff <- d.lincomb.dirichlet(x, a, alpha, method = method)
			# ans <- x^y * (1-x)^(m-y) * ff

			log.ff <- d.lincomb.dirichlet(x, a, alpha, method = method, log = TRUE)
			log.integral <- y*log(x) + (m-y)*log(1-x) + log.ff
			exp(log.integral)
		}

		fc <- numeric(J)

		for (j in 1:J)
		{
			a <- V[j,]

			tryCatch({
				II <- integrate(g, lower = min(a), upper = max(a), a = a,
					alpha = alpha, subdivisions = 400)
				fc[j] <- II$value
			},
			error = function(e)
			{
				warning("WARNING: R integrate function failed to converge. Trying trapezoidal rule...\n")
				II <- integrate.trap(g, xlim = sort(a), n = 100000, a = a, alpha = alpha)
				fc[j] <- II
			})
		}

		ff <- choose(m,y) * sum(Pi * fc)
	}
	else if (method == "mc")
	{
		lambda <- rdirichlet(mc.reps, alpha)
		mu <- V %*% t(lambda)

		fc <- mu^y * (1-mu)^(m-y)
		ff.mc <- choose(m,y) * t(fc) %*% Pi
		ff <- mean(ff.mc)
	}
	else
	{
		stop("Unknown method given in d.mixture.link.one")
	}

	if (log) return(log(ff))
	else return(ff)
}
d.mixture.link <- Vectorize(d.mixture.link.one, vectorize.args = c("y", "m", "p"))

# ------------------------- Likelihood -------------------------
likelihood.x <- function(y, X, m, Beta, Pi, kappa, log = FALSE, method = "imhof")
{
	p <- plogis(X %*% Beta)
	likelihood(y, m, p, Pi, kappa, log, method)
}

likelihood <- function(y, m, p, Pi, kappa, log = FALSE, method = "imhof")
{
	J <- length(Pi)
	n <- length(y)

	ff <- numeric(n)
	for (i in 1:n)
	{
		ff[i] <- d.mixture.link.one(y[i], m[i], p[i], Pi, kappa,
			log = log, method = method)
	}

	if (log) return(sum(ff))
	else return(prod(ff))
}

