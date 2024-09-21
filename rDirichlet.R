rDirichlet <-
function (n, alpha) 
{
    l <- length(alpha)
    x <- matrix(nrow = n, ncol = l)
    for (i in 1:l) {
        x[, i] <- rGamma(n, alpha[i])
    }
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}
