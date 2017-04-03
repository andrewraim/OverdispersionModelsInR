library(combinat)

dm.scoring <- function(k)
{
	theta.vec2list <- function(theta)
	{
		idx <- seq(1, k-1)
		Pi <- numeric(k)
		Pi[-k] <- theta[idx]
		Pi[k] <- 1 - sum(Pi)

		list(Pi = Pi, rho = theta[k])
	}
	
	theta.list2vec <- function(theta)
	{
		c(theta$Pi[-k], theta$rho)
	}

	loglik <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		sum(d.dm(Data$x, par$Pi, par$rho, Data$m, log = TRUE))
	}

	score <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		rho <- par$rho
		n <- Data$n
		m <- Data$m
		x <- Data$x
		if (length(m) == 1) m <- rep(m,n)

		C <- rho^(-2) * (1 - rho^2)
		alpha <- C * Pi
		h <- digamma

		score.part <- matrix(NA, n, k)
		for (j in 1:(k-1))
		{
			score.part[,j] <- C *
				(h(x[j,] + C*Pi[j]) - h(C*Pi[j]) - h(x[k,] + C*Pi[k]) + h(C*Pi[k]))
		}

		Pi.ext <- matrix(Pi, k, n)
		A <- h(x + C*Pi.ext) - h(C*Pi.ext)
		score.part[,k] <- -2*rho^(-3) * (h(C) - h(C+m) + t(A) %*% Pi)

		return(colSums(score.part))
	}

	approx.fim <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		rho <- par$rho
		M <- sum(Data$m)
		n <- Data$n

		C <- rho^(-2) * (1 - rho^2)
		alpha <- C * Pi
		afim <- matrix(NA, k, k)

		## Diagonal
		idx <- cbind(1:(k-1), 1:(k-1))
		afim[idx] <- n * C^2 * (trigamma(C*Pi[-k]) + trigamma(C*Pi[k]))

		## Lower right diagonal element
		afim[k,k] <- n * (sum(Pi^2 * trigamma(C*Pi)) - trigamma(C)) * 4 / rho^6

		## Off-diagonal
		for (i in 1:(k-1))
			for (j in seq(1, i-1, length.out = i-1)) {
				afim[i,j] <- n * C^2 * trigamma(C*Pi[k])
				afim[j,i] <- afim[i,j]
			}
		
		## Last row and column, except the bottom-right element
		idx <- 1:(k-1)
		afim[idx,k] <- n * C * (Pi[k]*trigamma(C*Pi[k]) - Pi[-k]*trigamma(C*Pi[-k])) * 2/rho^3
		afim[k,idx] <- afim[idx,k]

		return(afim)
	}

	fim <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		rho <- par$rho
		m <- Data$m
		n <- Data$n

		tab <- table(m)
		m.levels <- as.integer(names(tab))
		m.freqs <- as.integer(tab)

		qq <- k
		E <- matrix(0, qq, qq)
		
		for (i in 1:length(m.levels))
		{
			E.i <- matrix(0, qq, qq)
			out <- xsimplex(k, m.levels[i])
			
			for (j in 1:ncol(out))
			{
				z <- as.matrix(out[,j])
				u <- score(theta, Data = list(x = z, m = m.levels[i], n = 1))
				E.i <- E.i + u %*% t(u) * d.dm(z, Pi, rho, m.levels[i])
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
		
		qq <- k
		E <- matrix(0, qq, qq)
		for (r in 1:reps)
		{
			x <- r.dm(n, par$Pi, par$rho, m)
			S <- score(theta, Data = list(x = x, m = m, n = n))
			E <- E + S %*% t(S)
		}
		
		return(E / reps)
	}
	
	description <- "X[i,] ~ind~ DM_k(m[i], Pi, rho)"
	ret <- list(scoring.family = "dm", score = score, fim = fim,
		approx.fim = approx.fim, description = description, loglik = loglik,
		theta.vec2list = theta.vec2list, theta.list2vec = theta.list2vec,
		fim.mc = fim.mc)
	structure(ret, class = "scoring.family")
}
