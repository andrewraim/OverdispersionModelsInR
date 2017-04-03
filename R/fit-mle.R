fit.mle <- function(phi.init, loglik, theta.tx, extra.tx, Data, psi.names = NULL)
{
	stopifnot(!is.null(Data$n))
	n <- Data$n

	# Combine the two lists: theta.tx(phi) and extra.tx(phi)
	# into a vector. If the user provided psi.names, use those for the
	# variable names.
	psi.tx <- function(phi)
	{
		theta <- theta.tx(phi)
		psi1 <- unlist(theta)
		psi2 <- unlist(extra.tx(theta))
		psi <- c(psi1, psi2)
		
		if (!is.null(psi.names)) names(psi) <- psi.names
		return(psi)
	}

	optim.res <- optim(par = phi.init, fn = loglik, method = "L-BFGS-B",
		control = list(fnscale = -1, trace = 0), hessian = TRUE, Data = Data)

	phi.hat <- optim.res$par
	theta.hat <- theta.tx(phi.hat)
	xi.hat <- extra.tx(theta.hat)
	psi.hat <- psi.tx(phi.hat)

	V.phi <- -solve(optim.res$hessian)
	J.tx <- jacobian(psi.tx, phi.hat)
	V.psi <- J.tx %*% V.phi %*% t(J.tx)
	rownames(V.psi) <- colnames(V.psi) <- names(psi.hat)

	loglik.hat <- optim.res$value
	qq <- length(unlist(theta.hat))
	aic <- -2 * loglik.hat + 2*qq
	aicc <--2 * loglik.hat + 2*qq*n / (n-qq-1) 
	bic <- -2 * loglik.hat + qq*log(n)

	# Note: follow SAS NLMIXED, treat Est/SE as t with df = n
	df <- n
	se <- sqrt(diag(V.psi))
	t.val <- psi.hat / se
	p.val <- 2 * (1 - pt(abs(t.val), df = df))
	gr <- J.tx %*% grad(loglik, x = phi.hat, Data = Data)

	estimates <- cbind(psi.hat, se, t.val, p.val, gr)
	colnames(estimates) <- c("Estimate", "SE", "t-val", "P(|t|>t-val)", "Gradient")

	res <- list(estimates = estimates, loglik = loglik.hat, aic = aic,
		aicc = aicc, bic = bic, vcov = V.psi, optim.res = optim.res, df = df,
		qq = qq, description = "<Default>", theta.hat = theta.hat,
		xi.hat = xi.hat)
	class(res) <- "mle.fit"
	return(res)
}

confint.mle.fit <- function(fit.out, level = 0.95)
{
	dim.theta <- length(unlist(fit.out$theta.hat))
	dim.xi <- length(unlist(fit.out$xi.hat))

	na <- rep(NA, dim.theta + dim.xi)
	DF <- data.frame(Estimate = na, SE = na, Lower = na, Upper = na)
	rownames(DF) <- rownames(fit.out$estimates)

	w <- -qt((1-level)/2, df = fit.out$df)
	psi.hat <- fit.out$estimates[,1]
	se.psi.hat <- fit.out$estimates[,2]
	DF$Lower <- psi.hat - w * se.psi.hat
	DF$Upper <- psi.hat + w * se.psi.hat
	DF$Estimate <- psi.hat
	DF$SE <- se.psi.hat

	res <- list(ci = DF, df = fit.out$df, t.quantile = w, fit = fit.out, level = level)
	class(res) <- "mle.fit.ci"
	return(res)
}

print.mle.fit.ci <- function(ci.out)
{
	fit.out <- ci.out$fit
	dim.theta <- length(unlist(fit.out$theta.hat))
	dim.xi <- length(unlist(fit.out$xi.hat))

	printf("--- Parameter CIs (level %f) ---\n", ci.out$level)
	idx <- 1:dim.theta
	print(ci.out$ci[idx,])

	if (dim.xi > 0)
	{
		printf("--- Additional CIs (level %f) ---\n", ci.out$level)
		idx <- 1:dim.xi + dim.theta
		print(ci.out$ci[idx,])
	}
	
	printf("--\n")
	printf("Degrees of freedom = %d\n", fit.out$df)
	printf("t-quantile = %f\n", ci.out$t.quantile)
}

print.mle.fit <- function(fit.out)
{
	printf("Fit for model:\n")
	printf("%s\n", fit.out$description)

	dim.theta <- length(unlist(fit.out$theta.hat))
	dim.xi <- length(unlist(fit.out$xi.hat))

	DF <- as.data.frame(fit.out$estimates)
	DF[,1] <- round(DF[,1], 4)
	DF[,2] <- round(DF[,2], 4)
	DF[,3] <- round(DF[,3], 4)
	DF[,4] <- my.numerical.format(DF[,4])
	DF[,5] <- my.numerical.format(DF[,5])

	printf("--- Parameter Estimates ---\n")
	idx <- 1:dim.theta
	print(DF[idx,])

	if (dim.xi > 0)
	{
		printf("--- Additional Estimates ---\n")
		idx <- 1:dim.xi + dim.theta
		print(DF[idx,])
	}

	printf("--\n")
	printf("Degrees of freedom = %d\n", fit.out$df)
	printf("LogLik = %0.4f\n", fit.out$loglik)
	printf("AIC = %0.4f\n", fit.out$aic)
	printf("AICC = %0.4f\n", fit.out$aicc)
	printf("BIC = %0.4f\n", fit.out$bic)
}

coef.mle.fit <- function(fit.out)
{
	fit.out$estimates[,"Estimate"]
}
