gof.count.obs <- function(y, gof.breaks)
{
	Ob <- table(cut(y, gof.breaks, include.lowest = TRUE, right = FALSE))
	return(Ob)
}

# f is a density of the form f(x,i), for sample point x and the ith observation
# qq is the number of estimated parameters in f
gof.count <- function(y, f, gof.breaks, qq, max.y = 10000)
{
	K <- length(gof.breaks) - 1
	Ob <- gof.count.obs(y, gof.breaks)
	Ex <- numeric(K)
	n <- length(y)

	# This is a workaround since, without it, findInterval doesn't always seem to
	# be honoring the rightmost.closed = TRUE option, and counts end up in the
	# wrong intervals
	# NOTE that this seems to not be needed on dv410
	# gof.breaks[-c(1,K+1)] <- gof.breaks[-c(1,K+1)] + .Machine$double.eps

	tt <- seq(0, max.y)

	for (i in 1:n)
	{
		if (i %% 1000 == 0) printf("Computing expected count for obs %d\n", i)

		# Note: labelling of res columns (Int and ff) seems to get messed up if
		# ff is a matrix rather than a vector, so we do a coercion here.
		ff <- as.numeric( f(tt, i) )
		Int <- findInterval(x = tt, vec = gof.breaks, rightmost.closed = TRUE)
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
		breaks = gof.breaks, qq = qq)
	class(res) <- "gof.count"
	return(res)
}

print.gof.count <- function(gof.out)
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

plot.gof.count <- function(gof.out, ...)
{
	K <- nrow(gof.out$tab)
	coords <- barplot(gof.out$tab[,"Ob"], names.arg = 1:K,
		xlab = "Interval", ylab = "Frequency", ...)
	points(coords, gof.out$tab[,"Ex"], pch = 19)
	box()
}

