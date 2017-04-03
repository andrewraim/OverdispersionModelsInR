gof.binomial.obs <- function(y, m, gof.breaks)
{
	Ob <- table(cut(y/m, gof.breaks, include.lowest = TRUE))
	return(Ob)
}

# f is a density of the form f(x,i), for sample point x and the ith observation
# qq is the number of estimated parameters in f
gof.binomial <- function(y, m, f, gof.breaks, qq)
{
	if (length(m) == 1) m <- rep(m, length(y))
	K <- length(gof.breaks) - 1
	Ob <- gof.binomial.obs(y, m, gof.breaks)
	Ex <- numeric(K)
	n <- length(y)

	# This is a workaround since, without it, findInterval doesn't always seem to
	# be honoring the rightmost.closed = TRUE option, and counts end up in the
	# wrong intervals
	gof.breaks[-c(1,K+1)] <- gof.breaks[-c(1,K+1)] + .Machine$double.eps

	for (i in 1:n)
	{
		# Note: labelling of res columns (Int and ff) seems to get messed up if
		# ff is a matrix rather than a vector, so we do a coercion here.
		tt <- seq(0, m[i])
		ff <- as.numeric( f(tt, i) )
		Int <- findInterval(x = tt / m[i], vec = gof.breaks, rightmost.closed = TRUE)
		res <- aggregate(ff ~ Int, FUN = sum)
		Ex[res$Int] <- Ex[res$Int] + res$ff
	}

	tab <- cbind(Ob, Ex)
	X <- sum( (Ob - Ex)^2 / Ex )
	df.low <- K - 1 - qq
	df.high <- K - 1
	pvalue.low <- 1 - pchisq(X, df.low)
	pvalue.high <- 1 - pchisq(X, df.high)
		
	brak <- c("[", rep("(", K-1))
	labels <- sprintf("%s%0.04f,%0.04f]", brak, gof.breaks[-(K+1)], gof.breaks[-1])
	rownames(tab) <- labels 
		
	res <- list(tab = tab, X = X, df.low = df.low, df.high = df.high,
		pvalue.low = pvalue.low, pvalue.high = pvalue.high, f = f,
		breaks = gof.breaks, m = m, qq = qq)
	class(res) <- "gof.binomial"
	return(res)
}

print.gof.binomial <- function(gof.out)
{
	printf("GOF for model:\n")
	print(gof.out$f)
	printf("--\n")
	print(round(gof.out$tab, 4))
	printf("--\n")
	printf("X = %0.4f\n", gof.out$X)
	printf("DF in [%d, %d]\n", gof.out$df.low, gof.out$df.high)
	printf("p-value in [%f, %f]\n", gof.out$pvalue.low, gof.out$pvalue.high)
}

plot.gof.binomial <- function(gof.out, ...)
{
	K <- nrow(gof.out$tab)
	coords <- barplot(gof.out$tab[,"Ob"], names.arg = 1:K,
		xlab = "Interval", ylab = "Frequency", ...)
	points(coords, gof.out$tab[,"Ex"], pch = 19)
	box()
}
