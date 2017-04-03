fit.mixture.link.mle <- function(y, m, J, extra.tx = null.tx,
	var.names = NULL, method = "imhof", phi.init = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))

	Data <- list(y = y, m = m, n = length(y))
	qq <- 2 + J
	
	if (is.null(phi.init))
	{
		phi.init <- c(0, qlogis(normalize(1:J)), log(1))
	}
	
	theta.tx <- function(phi)
	{
		Pi <- normalize(plogis(phi[1:J + 1]))
		list(p = plogis(phi[1]), Pi = Pi, kappa = exp(phi[qq]))
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		
		# cat("p = ", theta$p)
		# cat(", Pi = (", theta$Pi, ")")
		# cat(", kappa = ", theta$kappa)
		# cat("\n")
		
		ll <- numeric(Data$n)
		for (i in 1:Data$n)
		{
			ll[i] <- d.mixture.link.one(Data$y[i], Data$m[i], theta$p, theta$Pi,
				theta$kappa, method = method, log = TRUE)
		}

		return(sum(ll))
	}

	start <- Sys.time()
	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- paste0(
		sprintf("y[i] ~iid~ MixLink_%d(m[i], p, Pi, kappa),\n", J),
		sprintf("Using method: %s\n", method),
		sprintf("Elapsed time: %0.02f sec", as.numeric(Sys.time() - start, units = "secs"))
	)
	return(fit.out)
}

fit.mixture.link.x.mle <- function(y, m, X, J, extra.tx = null.tx,
	var.names = NULL, method = "imhof", phi.init = NULL)
{
	if (length(m) == 1) m <- rep(m, length(y))
	
	Data <- list(y = y, m = m, X = X, n = length(y))
	d <- ncol(X)
	qq <- d + J + 1

	if (is.null(phi.init))
	{
		phi.init <- c(rep(0,d), qlogis(normalize(1:J)), log(1))
	}

	theta.tx <- function(phi)
	{
		list(
			Beta = phi[1:d],
			Pi = normalize(plogis(phi[1:J + d])),
			kappa = exp(phi[qq])
		)
	}
	
	loglik <- function(phi, Data)
	{
		theta <- theta.tx(phi)
		
		cat("Beta = (", theta$Beta, ")")
		cat(", Pi = (", theta$Pi, ")")
		cat(", kappa = ", theta$kappa)
		cat("\n")
		
		p <- plogis(X %*% theta$Beta)
		ll <- numeric(Data$n)
		for (i in 1:Data$n)
		{
			ll[i] <- d.mixture.link.one(Data$y[i], Data$m[i], p[i], theta$Pi,
				theta$kappa, method = method, log = TRUE)
		}
		
		return(sum(ll))
	}

	start <- Sys.time()
	fit.out <- fit.mle(phi.init, loglik, theta.tx, extra.tx,
		Data, psi.names = var.names)
	
	fit.out$description <- fit.out$description <- paste0(
		sprintf("y[i] ~indep~ MixLink_%d(m[i], p[i], Pi, kappa)\n", J),
		"logit( E(Y[i]) ) = x[i]^T Beta\n",
		sprintf("Using method: %s\n", method),
		sprintf("Elapsed time: %f sec", as.numeric(Sys.time() - start, units = "secs"))
	)
	return(fit.out)
}
