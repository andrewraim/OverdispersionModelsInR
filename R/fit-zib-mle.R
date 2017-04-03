fit.zib.mle <- function(y, m, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))
	
	Data <- list(y = y, m = m, n = length(y))
	qq <- 2
	
	phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(Pi = plogis(phi[1]), rho = plogis(phi[2]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		sum( d.zib(Data$y, Pi = theta$Pi, rho = theta$rho,
			m = Data$m, log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- "y[i] ~iid~ ZIB(m[i], Pi, rho)"
	return(fit.out)
}

fit.zib.x.mle <- function(y, m, X, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, X = X, n = length(y))
	d <- ncol(X)
	qq <- d + 1
	
	phi.init <- rep(0, qq)

	theta.tx <- function(phi)
	{
		list(Beta = phi[1:d], rho = plogis(phi[d+1]))
	}

	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		Pi <- plogis(Data$X %*% theta$Beta)
		sum( d.zib(Data$y, Pi = Pi, rho = theta$rho, m = Data$m, log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- paste0(
		"y[i] ~indep~ ZIB(m[i], Pi[i], rho)\n",
		"logit(Pi[i]) = x[i]^T Beta")
	return(fit.out)
}

fit.zib.xz.mle <- function(y, m, X, Z, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))
	
	Data <- list(y = y, m = m, X = X, Z = Z, n = length(y))
	d1 <- ncol(X)
	d2 <- ncol(Z)
	qq <- d1 + d2	
	
	phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(Beta = phi[1:d1], Gamma = phi[1:d2 + d1])
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		Pi <- plogis(Data$X %*% theta$Beta)
		rho <- plogis(Data$Z %*% theta$Gamma)
		sum( d.rcb(Data$y, Pi = Pi, rho = rho, m = Data$m, log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		"y[i] ~indep~ ZIB(m[i], Pi[i], rho[i])\n",
		"logit(Pi[i]) = x[i]^T Beta\n",
		"logit(rho[i]) = z[i]^T Gamma")
	return(fit.out)
}

