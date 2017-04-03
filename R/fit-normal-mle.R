fit.normal.mle <- function(y, extra.tx = null.tx, var.names = NULL)
{
	Data <- list(y = y, n = length(y))
	qq <- 2
	
	phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(mu = phi[1], sigma = exp(phi[2]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		sum( dnorm(Data$y, mean = theta$mu, sd = theta$sigma, log = TRUE) )
	}
	
	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- "y[i] ~iid~ N(mu, sigma^2)"
	return(fit.out)
}

fit.normal.x.mle <- function(y, X, extra.tx = null.tx, var.names = NULL)
{
	Data <- list(y = y, X = X, n = length(y))
	d <- ncol(X)
	qq <- d + 1
	
	phi.init <- rep(0, qq)
	
	theta.tx <- function(phi)
	{
		list(Beta = phi[1:d], sigma = exp(phi[d+1]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		mu <- Data$X %*% theta$Beta
		sum( dnorm(Data$y, mean = mu, sd = theta$sigma, log = TRUE) )
	}
	
	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		"y[i] ~indep~ N(mu[i], sigma^2)\n",
		"mu[i] = x[i]^T Beta")
	return(fit.out)
}
