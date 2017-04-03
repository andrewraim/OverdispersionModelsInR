fit.dm.mle <- function(y, m, extra.tx = null.tx, var.names = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))
	k <- nrow(y)

	Data <- list(y = y, m = m, n = length(y))
	qq <- k

	phi.init <- rep(0, qq)

	theta.tx <- function(phi)
	{
		idx <- 1:(k-1)
		list(Pi = prob.simplex.tx(phi[idx]), rho = plogis(phi[k]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		sum( d.dm(Data$y, Pi = theta$Pi, rho = theta$rho,
			m = Data$m, log = TRUE) )
	}

	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- sprintf("y[i] ~iid~ DM_%d(m[i], Pi, rho)", k)
	return(fit.out)
}
