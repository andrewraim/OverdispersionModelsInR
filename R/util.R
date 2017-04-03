printf <- function(msg, ...)
{
	cat(sprintf(msg, ...))
}

fprintf <- function(file, msg, ...)
{
	cat(sprintf(msg, ...), file = file)
}

my.numerical.format <- function(x, lower = 1e-4)
{
	idx1 <- which(abs(x) < lower)
	idx2 <- setdiff(1:length(x), idx1)
	y <- character(length(x))
	y[idx1] <- sprintf("%0.3E", x[idx1])
	y[idx2] <- sprintf("%0.4f", x[idx2])
	names(y) <- names(x)
	return(y)
}

print.debug <- function(x, debug = FALSE)
{
	if (debug) print(x)
}

normalize <- function(x) { x / sum(x) }

null.tx <- function(phi) { list() }

# Transform from probability simplex S^J to R^(J-1)
mlogit <- function(p)
{
        J <- length(p)
        x <- log(p[-J] / p[J])
        return(x)
}

# Transform from R^(J-1) to probability simplex S^J
inv.mlogit <- function(x)
{
        z <- exp(x)
        P.J <- 1 / (1 + sum(z))
        p <- c(z * P.J, P.J)
        return(p)
}

integrate.trap <- function(f, xlim, n, ...)
{
	a <- xlim[1]
	b <- xlim[2]
	
	i <- seq(1, n)
	ff <- f(a + i/2 * (b-a)/n, ...)
	
	II <- (b - a) / (2*n) * sum(ff)
	return(II)
}


rpad.str <- function(x, len)
{
	L <- length(x)
	y <- character(L)
	pad <- len - nchar(x)

	for (i in 1:L)
	{
		s <- paste0(rep(" ", pad[i]), collapse="")
		y[i] <- sprintf("%s%s", x[i], s)
	}

	return(y)
}


## The following numerical derivative functions are not currently needed.
## We'll use numDeriv functions instead.

if (FALSE)
grad <- function(func, x,  method.args = list(d = 1e-04), ...)
{
	k <- length(x)
	eps <- method.args$d
	E <- diag(eps, k) 
  
	fx <- func(x, ...)
	D <- numeric(k)
  
	for (j in 1:k)
	{
		D[j] <- func(x + E[j,], ...) - fx
	}
  
	return(D / eps)
}

if (FALSE)
hessian <- function(func, x, method.args = list(d = 1e-04), ...)
{
	k <- length(x)
	eps <- method.args$d
	E <- diag(eps, k) ## Perturbations in each axis direction
  
	fx <- func(x, ...)
	H <- matrix(NA, k, k)
  
	for (i in 1:k)
	{
		for (j in 1:k)
		{ 
			H[i,j] <- func(x + E[i,] + E[j,], ...) - 
				func(x + E[i,], ...) -
				func(x + E[j,], ...) + fx
		}
	}
  
	return(H / eps^2)
}

if (FALSE)
jacobian <- function(func, x, method.args = list(d = 1e-04), ...)
{
	k <- length(x)
	eps <- method.args$d
	E <- diag(eps, k) ## Perturbation in each axis direction
  
	fx <- func(x, ...)
	n <- length(fx)
	J <- matrix(NA, n, k)
  
	for (j in 1:k)  
	{
		fxe <- func(x + E[j,], ...)
		J[,j] <- fxe - fx
	}
  
	return(J / eps)
}
