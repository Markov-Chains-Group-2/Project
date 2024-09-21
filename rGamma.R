rGamma <-
function (N, SHAPE, RATE = 1, SCALE = 1/RATE) 
{
    RGAMMA <- function(n, Shape, rate = 1, scale = 1/rate) {
        scale1 <- function(shape) {
            t <- 0.07 + 0.75 * sqrt(1 - shape)
            b <- 1 + (shape * exp(-t))/t
            YY <- 0
            while (YY == 0) {
                u1 <- runif(1)
                u2 <- runif(1)
                v <- b * u1
                if (v <= 1) {
                  x <- t * (v^(1/shape))
                  if (u2 <= (2 - x)/(2 + x)) {
                    YY <- x
                  }
                  else {
                    if (u2 <= exp(-x)) {
                      YY <- x
                    }
                  }
                }
                else {
                  x <- -log(t * (b - v)/shape)
                  y <- x/t
                  if (u2 * (shape + y * (1 - shape)) <= 1) {
                    YY <- x
                  }
                  else {
                    if (u2 <= y^(shape - 1)) {
                      YY <- x
                    }
                  }
                }
            }
            return(YY)
        }
        dummy <- 0
        for (i in 1:n) {
            dummy2 <- scale1(Shape)
            dummy[i] <- dummy2 * scale
        }
        return(dummy)
    }
    if (SHAPE <= 1) {
        return(RGAMMA(n = N, Shape = SHAPE, scale = SCALE))
    }
    else (return(rgamma(N, shape = SHAPE, scale = SCALE)))
}
