fit.binmix.mle <- function(y, m, J, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, n = length(y))
	qq <- 2*J

	phi.init <- phi.init <- c((1:J)/(J+1), normalize(1:J))

	theta.tx <- function(phi)
	{
		idx <- 1:J + J
		Pi <- normalize(plogis(phi[idx]))
		list(p = plogis(phi[1:J]), Pi = Pi)
	}

	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		sum(d.binom.mix(Data$y, theta$p, theta$Pi, Data$m, log = TRUE))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- sprintf("y[i] ~indep~ BinMix_%d(m[i], p, Pi)", J)
	return(fit.out)
}

## X.g should be a list of covariate matrices of length J
fit.binmix.x.mle <- function(y, m, X.g, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, X.g = X.g, n = length(y))
	d <- unlist(Map(ncol, X.g))
	J <- length(d)
	qq <- sum(d) + J

	phi.init <- c(rep(0, sum(d)), normalize(1:J))

	theta.tx <- function(phi)
	{
		idx <- 1:sum(d)
		grp <- as.integer( mapply(rep, times = d, x = 1:length(d)) )
		Beta.g <- split(phi[idx], grp)

		idx <- seq(1,J) + sum(d)
		Pi <- normalize(plogis(phi[idx]))

		list(Beta.g = Beta.g, Pi = Pi)
	}

	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		ff <- matrix(NA, Data$n, J)

		for (j in 1:J)
		{
			p.j <- plogis(Data$X.g[[j]] %*% theta$Beta.g[[j]])
			ff[,j] <- dbinom(Data$y, p = p.j, size = Data$m)
		}
		
		sum(log(ff %*% theta$Pi))
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- paste0(
		sprintf("y[i] ~indep~ BinMix_%d(m[i], p[i], Pi)\n", J),
		"logit(p[i,j]) = x[i]^T Beta[j]")
	return(fit.out)
}
