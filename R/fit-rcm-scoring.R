library(combinat)

rcm.scoring <- function(k)
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
		sum(d.rcm(Data$x, par$Pi, par$rho, Data$m, log = TRUE))
	}

	score <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		rho <- par$rho
		m <- Data$m
		n <- Data$n

		## Compute the score of MultMix_k with approproiate parameters, then
		## transform to RCM score

		II <- diag(k-1)

		eta <- numeric((k-1)*(k+1))
		for (j in 1:(k-1))
		{
			idx <- 1:(k-1) + (k-1)*(j-1)
			eta[idx] <- (1-rho) * Pi[-k] + rho*II[,j]
		}
		
		idx <- 1:(k-1) + (k-1)*(k-1)
		eta[idx] <- (1-rho) * Pi[-k]

		idx <- 1:(k-1) + k*(k-1)
		eta[idx] <- Pi[-k]

		multmix <- multmix.scoring(k, k)
		S <- multmix$score(eta, Data)

		J <- matrix(NA, k, (k+1)*(k-1))
		for (j in 1:k)
		{
			idx1 <- 1:(k-1)
			idx2 <- 1:(k-1) + (j-1)*(k-1)
			J[idx1, idx2] <- (1-rho)*II
			
			if (j < k)
				J[k, idx2] <- -t(Pi[-k]) + t(II[,j])
			else
				J[k, idx2] <- -t(Pi[-k])
		}
		
		idx1 <- 1:(k-1)
		idx2 <- 1:(k-1) + k*(k-1)
		J[idx1, idx2] <- II

		J[k, idx2] <- 0

		return(as.numeric(J %*% S))
	}

	approx.fim <- function(theta, Data)
	{
		par <- theta.vec2list(theta)
		Pi <- par$Pi
		rho <- par$rho
		m <- Data$m
		n <- Data$n

		Beta <- Pi / ((1 - rho)*Pi + rho) + (1-Pi) / ((1 - rho)*Pi)
		Gamma <- Pi*(1 - Pi) / ((1 - rho)*Pi + rho) + Pi / (1 - rho)
		M <- sum(m)
		
		afim <- matrix(NA, k, k)
		
		## Diagonal
		for (i in 1:(k-1)) {
			afim[i,i] <- M * (1 - rho)^2 * (Beta[i] + Beta[k]) - (1/Pi[i] + 1/Pi[k])
		}

		## Lower right diagonal element
		afim[k,k] <- M/(1 - rho) * sum(Pi*(1-Pi) / ((1 - rho)*Pi + rho) )
		
		## Off-diagonal
		for (i in 1:(k-1))
		for (j in seq(1, i-1, length.out = i-1)) {
			afim[i,j] <- M * (1 - rho)^2 * Beta[k] - 1/Pi[k]
			afim[j,i] <- afim[i,j]
		}

		## Last row and column, except the bottom-right element
		for (i in 1:(k-1)) {
			afim[i,k] <- M * (1 - rho) * (Gamma[i] - Gamma[k])
			afim[k,i] <- afim[i,k]
		}

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
				E.i <- E.i + u %*% t(u) * d.rcm(z, Pi, rho, m.levels[i])
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
			x <- r.rcm(n, par$Pi, par$rho, m)
			S <- score(theta, Data = list(x = x, m = m, n = n))
			E <- E + S %*% t(S)
		}
		
		return(E / reps)
	}
	
	description <- "X[i,] ~ind~ RCM_k(m[i], Pi, rho)"
	ret <- list(scoring.family = "rcm", score = score, fim = fim,
		approx.fim = approx.fim, description = description, loglik = loglik,
		theta.vec2list = theta.vec2list, theta.list2vec = theta.list2vec,
		fim.mc = fim.mc)
	structure(ret, class = "scoring.family")
}
