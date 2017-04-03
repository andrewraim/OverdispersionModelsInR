library(combinat)

multmix.scoring <- function(s, k)
{
	theta.vec2list <- function(theta)
	{
		idx <- seq(1, s*(k-1))
		P <- matrix(0, k, s)
		P[-k,] <- matrix(theta[idx], k-1, s, byrow = FALSE)
		P[k,] <- 1 - colSums(P)

		idx <- seq(1, s-1) + s*(k-1)
		Pi <- numeric(s)
		Pi[-s] <- theta[idx]
		Pi[s] <- 1 - sum(Pi)

		list(P = P, Pi = Pi)
	}

	theta.list2vec <- function(theta)
	{
		c(as.numeric(theta$P[-k,]), theta$Pi[-s])
	}

	loglik <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		sum(d.multmix(Data$x, par$P, par$Pi, Data$m, log = TRUE))
	}

	score <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		P <- par$P
		x <- as.matrix(Data$x)
		m <- Data$m
		n <- Data$n
		qq <- s*k - 1

		U <- matrix(NA, n, qq)
		for (i in 1:n)
		{
			w <- d.multmix(x[,i], P, Pi, m[i])
			u <- numeric(qq)

			for (l in 1:s)
			{
				Z <- x[-k,i]/P[-k,l] - x[k,i]/P[k,l]
				v <- dmultinom(x[,i], prob = P[,l], size = m[i])
				u_l <- Pi[l] * Z * v/w

				a <- (l-1) * (k-1) + 1
				idx <- seq(a, a+k-2)
				u[idx] <- u_l
			}

			Y <- numeric(s-1)
			for (l in 1:(s-1))
			{
				Y[l] <- dmultinom(x[,i], prob = P[,l], size = m[i]) -
					dmultinom(x[,i], prob = P[,s], size = m[i])
			}
			u_Pi <- Y / w
			idx <- seq( (k-1)*s + 1, qq )
			u[idx] <- u_Pi

			U[i,] <- u
		}

		return(colSums(U))
	}

	approx.fim <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- as.matrix(par$Pi)
		P <- par$P
		m <- Data$m
		n <- Data$n

		qq <- s*k-1
		fim <- matrix(0, qq, qq)
		ones <- as.matrix(rep(1, k-1))

		for (l in 1:s)
		{
			a <- (l-1) * (k-1) + 1
			idx <- seq(a, a+k-2)
			D.inv <- diag(1/P[-k,l], k-1)
			fim[idx,idx] <- Pi[l] * sum(m) * (D.inv + 1/P[k,l] * ones %*% t(ones))
		}

		idx <- seq( (k-1)*s + 1, qq )
		ones <- as.matrix(rep(1, s-1))
		D.inv <- diag(1/Pi[-s], s-1)
		fim[idx,idx] <- n * (D.inv + 1/Pi[s] * ones %*% t(ones))

		return(fim)
	}

	fim <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- as.matrix(par$Pi)
		P <- par$P
		m <- Data$m
		n <- Data$n

		tab <- table(m)
		m.levels <- as.integer(names(tab))
		m.freqs <- as.integer(tab)

		# s <- length(Pi)
		# k <- nrow(P)
		qq <- s*k-1

		E <- matrix(0, qq, qq)

		for (i in 1:length(m.levels))
		{
			E.i <- matrix(0, qq, qq)
			out <- xsimplex(k, m.levels[i])

			for (j in 1:ncol(out))
			{
				z <- as.matrix(out[,j])
				u <- score(theta, Data = list(x = z, m = m.levels[i], n = 1))
				E.i <- E.i + u %*% t(u) * d.multmix(z, P, Pi, m.levels[i])
			}

			E <- E + m.freqs[i] * E.i
		}

		return(E)
	}

	fim.mc <- function(theta, Data, reps = 100000)
	{
		par <- theta.vec2list(theta)
		m <- Data$m
		n <- Data$n

		qq <- s*k-1
		E <- matrix(0, qq, qq)
		for (r in 1:reps)
		{
			y <- r.multmix(n, par$P, par$Pi, m)
			S <- score(theta, Data = list(x = y, m = m, n = n))
			E <- E + S %*% t(S)
		}

		return(E / reps)
	}
	
	description <- "X[i,] ~ind~ MultMix_s(m[i], P, Pi)"
	ret <- list(scoring.family = "multmix", score = score, fim = fim,
		approx.fim = approx.fim, description = description, loglik = loglik,
		theta.vec2list = theta.vec2list, theta.list2vec = theta.list2vec,
		fim.mc = fim.mc)
	structure(ret, class = "scoring.family")
}
