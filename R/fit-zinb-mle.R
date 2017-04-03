fit.zinb.mle <- function(y, extra.tx = null.tx, var.names = NULL,
	phi.init = NULL)
{
	Data <- list(y = y, n = length(y))
	qq <- 3
	
	if (is.null(phi.init))
		phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(mu = exp(phi[1]), Pi = plogis(phi[2]), kappa = exp(phi[3]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		sum( d.zinb(Data$y, mu = theta$mu, kappa = theta$kappa, Pi = theta$Pi,
			log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- "y[i] ~iid~ ZINB(mu, kappa, Pi)"
	return(fit.out)
}

fit.zinb.x.mle <- function(y, X, extra.tx = null.tx, var.names = NULL,
	phi.init = NULL)
{
	Data <- list(y = y, X = X, n = length(y))
	d <- ncol(X)
	qq <- d + 2
	
	if (is.null(phi.init))
		phi.init <- rep(0, qq)

	theta.tx <- function(phi)
	{
		list(Beta = phi[1:d], Pi = plogis(phi[d+1]), kappa = exp(phi[d+2]))
	}

	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		mu <- exp(Data$X %*% theta$Beta)
		sum( d.zinb(Data$y, mu = mu, kappa = theta$kappa, Pi = theta$Pi,
			log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)

	fit.out$description <- paste0(
		"y[i] ~indep~ ZINB(mu[i], kappa, Pi)\n",
		"exp(mu[i]) = x[i]^T Beta")
	return(fit.out)
}

fit.zinb.xz.mle <- function(y, X, Z, extra.tx = null.tx, var.names = NULL,
	phi.init = NULL)
{
	Data <- list(y = y, X = X, Z = Z, n = length(y))
	d1 <- ncol(X)
	d2 <- ncol(Z)
	qq <- d1 + d2 + 1
	
	if (is.null(phi.init))
		phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(Beta = phi[1:d1], Gamma = phi[1:d2 + d1], kappa = exp(phi[d1+d2+1]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		mu <- exp(Data$X %*% theta$Beta)
		Pi <- plogis(Data$Z %*% theta$Gamma)
		sum( d.zinb(Data$y, mu = mu, kappa = theta$kappa, Pi = Pi, log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		"y[i] ~indep~ ZINB(mu[i], kappa, Pi[i])\n",
		"exp(mu[i]) = x[i]^T Beta\n",
		"logit(Pi[i]) = z[i]^T Gamma")
	return(fit.out)
}
