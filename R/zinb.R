r.zinb <- function(n, mu, kappa, Pi)
{
    z <- rbinom(n, 1, Pi)
    x <- rnbinom(n, mu = mu, size = 1/kappa)
    (z==1)*0 + (z==0)*x
}

d.zinb <- function(x, mu, kappa, Pi, log = FALSE)
{
    fz <- (x == 0)
    f1 <- dnbinom(x, mu = mu, size = 1/kappa, log = FALSE)
    ff <- Pi*fz + (1-Pi)*f1
    if (log) return(log(ff))
    else return(ff)
}

