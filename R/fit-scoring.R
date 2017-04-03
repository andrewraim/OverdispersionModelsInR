# -------------------- "Raw" versions of NR-type algorithms --------------------
## User supplies all the arguments themselves

fit.nr <- function(theta.init, loglik, score = NULL, hess = NULL,
	Data, max.iter = Inf, tol = 1e-6, extra.tx = null.tx, var.names = NULL,
	alg.name = "Newton-Raphson", description = "<Default>")
{
	n <- Data$n
	ll <- -Inf
	iter <- 0
	delta <- Inf
	theta <- theta.init

	if (is.null(score))
	{
		printf("No gradient function was specified. Using a numerical one.\n")
		score <- function(theta, Data) { grad(func = loglik, x = theta, Data = Data,
			method.args = list(d = 1e-4)) }
	}

	if (is.null(hess))
	{
		printf("No hessian function was specified. Using a numerical one.\n")
		hess <- function(theta, Data) { hessian(func = loglik, x = theta, Data = Data,
			method.args = list(d = 1e-4)) }
	}

	while (iter < max.iter && abs(delta) > tol)
	{
		iter <- iter + 1
		S <- score(theta, Data)
		H <- hess(theta, Data)
		theta <- theta - solve(H, S)

		ll.old <- ll
		ll <- loglik(theta, Data)
		delta <- ll - ll.old
		printf("%s: After iter %d: loglik = %f, delta = %e, estimates:\n",
			alg.name, iter, ll, delta)
		print(theta)
	}

	# Combine the two lists: theta and extra.tx(theta) into a vector. If the
	# user provided var.names, use those for the variable names.
	psi.tx <- function(theta)
	{
		psi1 <- theta
		psi2 <- unlist(extra.tx(theta))
		psi <- c(psi1, psi2)
		
		if (!is.null(var.names)) names(psi) <- var.names
		return(psi)
	}

	theta.hat <- theta
	xi.hat <- extra.tx(theta.hat)
	psi.hat <- psi.tx(theta.hat)

	V.theta <- -solve(H)
	J.tx <- jacobian(psi.tx, theta.hat)
	V.psi <- J.tx %*% V.theta %*% t(J.tx)
	rownames(V.psi) <- colnames(V.psi) <- names(psi.hat)

	loglik.hat <- ll
	qq <- length(theta.hat)
	aic <- -2 * loglik.hat + 2*qq
	aicc <--2 * loglik.hat + 2*qq*n / (n-qq-1) 
	bic <- -2 * loglik.hat + qq*log(n)

	# Note: follow SAS NLMIXED, treat Est/SE as t with df = n
	df <- Data$n
	se <- sqrt(diag(V.psi))
	t.val <- psi.hat / se
	p.val <- 2 * (1 - pt(abs(t.val), df = df))
	gr <- J.tx %*% score(theta.hat, Data)

	estimates <- cbind(psi.hat, se, t.val, p.val, gr)
	colnames(estimates) <- c("Estimate", "SE", "t-val", "P(|t|>t-val)", "Gradient")

	res <- list(estimates = estimates, loglik = loglik.hat, aic = aic,
		aicc = aicc, bic = bic, vcov = V.psi, iter = iter, tol = delta,
		converged = (iter < max.iter), df = df, qq = qq,
		description = description, theta.hat = theta.hat,
		xi.hat = xi.hat)

	class(res) <- "nr.fit"
	return(res)
}

fit.fs <- function(theta.init, loglik, score, fim, Data, max.iter = Inf,
	tol = 1e-6, extra.tx = null.tx, var.names = NULL, description = "<Default>")
{
	hess <- function(theta, Data) { -fim(theta, Data) }
	fit.nr(theta.init, loglik, score, hess, Data, max.iter, tol,
		extra.tx, var.names, alg.name = "Fisher Scoring", description)
}

fit.afs <- function(theta.init, loglik, score, afim, Data, max.iter = Inf,
	tol = 1e-6, extra.tx = null.tx, var.names = NULL, description = "<Default>")
{
	hess <- function(theta, Data) { -afim(theta, Data) }
	fit.nr(theta.init, loglik, score, hess, Data, max.iter, tol,
		extra.tx, var.names, alg.name = "Approx Scoring", description)
}

fit.fs.hybrid <- function(theta.init, loglik, score, fim, afim, Data,
	max.iter = Inf, tol = 1e-6, warmup.tol = 1e-4, extra.tx = null.tx, var.names = NULL,
	description = "<Default>")
{
	printf("FS Hybrid: Doing Approx Scoring until warmup.tol = %g\n", warmup.tol)
	out.afs <- fit.afs(theta.init, loglik, score, afim, Data, max.iter,
		tol = warmup.tol)

	printf("FS Hybrid: Doing Fisher Scoring until tol = %g\n", tol)
	res <- fit.fs(out.afs$theta.hat, loglik, score, fim, Data, max.iter - out.afs$iter,
		tol, extra.tx, var.names, description)

	res$iter <- res$iter + out.afs$iter
	return(res)
}

fit.nr.hybrid <- function(theta.init, loglik, score, hess, afim, Data,
	max.iter = Inf, tol = 1e-6, warmup.tol = 1e-4, extra.tx = null.tx, var.names = NULL,
	description = "<Default>")
{
	printf("NR Hybrid: Doing Approx Scoring until warmup.tol = %g\n", warmup.tol)
	out.afs <- fit.afs(theta.init, loglik, score, afim, Data, max.iter, tol = warmup.tol)

	printf("NR Hybrid: Doing Newton-Raphson until tol = %g\n", tol)
	res <- fit.nr(out.afs$theta.hat, loglik, score, hess, Data, max.iter - out.afs$iter,
		tol, extra.tx, var.names, description = description)

	res$iter <- res$iter + out.afs$iter
	return(res)
}

# -------------------- "Family" versions of NR-type algorithms --------------------
## User specifies a scoring.family, which supplies all the functions needed for
## NR-type algorithms

fit.family.nr <- function(theta.init, family, Data, max.iter = Inf, tol = 1e-6,
					   extra.tx = null.tx, var.names = NULL, alg.name = "Newton-Raphson")
{
	stopifnot(class(family) == "scoring.family")
	fit.nr(theta.init, family$loglik, family$score, family$hess,
	   Data, max.iter, tol, extra.tx, var.names, alg.name, description = family$description)
}

fit.family.fs <- function(theta.init, family, Data, max.iter = Inf, tol = 1e-6,
	extra.tx = null.tx, var.names = NULL)
{
	stopifnot(class(family) == "scoring.family")
	fit.fs(theta.init, family$loglik, family$score, family$fim, Data,
		max.iter, tol, extra.tx, var.names, description = family$description)
}

fit.family.afs <- function(theta.init, family, Data, max.iter = Inf, tol = 1e-6,
	extra.tx = null.tx, var.names = NULL)
{
	stopifnot(class(family) == "scoring.family")
	fit.afs(theta.init, family$loglik, family$score, family$approx.fim, Data,
		   max.iter, tol, extra.tx, var.names, description = family$description)
}

fit.family.fs.hybrid <- function(theta.init, family, Data, max.iter = Inf, tol = 1e-6,
	warmup.tol = 1e-4, extra.tx = null.tx, var.names = NULL)
{
	stopifnot(class(family) == "scoring.family")
	fit.fs.hybrid(theta.init, family$loglik, family$score, family$fim,
		family$approx.fim, Data, max.iter, tol, warmup.tol, extra.tx, var.names,
		description = family$description)
}

fit.family.nr.hybrid <- function(theta.init, family, Data, max.iter = Inf, tol = 1e-6,
	warmup.tol = 1e-4, extra.tx = null.tx, var.names = NULL)
{
	stopifnot(class(family) == "scoring.family")
	fit.nr.hybrid(theta.init, family$loglik, family$score, family$hess,
		family$approx.fim, Data, max.iter, tol, warmup.tol, extra.tx, var.names,
		description = family$description)
}

# -------------------- Utility functions for NR-type algorithms --------------------
print.nr.fit <- function(fit.out)
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
	printf("Iterations = %d\n", fit.out$iter)
	printf("Tolerance = %g\n", fit.out$tol)
}

confint.nr.fit <- function(fit.out, level = 0.95)
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
	class(res) <- "nr.fit.ci"
	return(res)
}
