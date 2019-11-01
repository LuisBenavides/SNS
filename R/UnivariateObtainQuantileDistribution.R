#' @title Obtain Quantile from Distribution Function
#' @description Get the quantile \code{theta} from several distributions with user defined mean and variance.
#' @inheritParams getDist
#' @inheritParams NS
#' @return A quantile \code{theta} of the selected \code{Ftheta} distribution with its parameters.
#' @export
#' @examples
#' getQuantile("Normal", c(0,1), 0, 1, 0.5)
getQuantile <- function(mu, sigma, Ftheta, dist,
                        par.location = 0, par.scale = 1, par.shape = 1, dist.par = NULL) {
  switch(dist,
         Uniform  = {
           a <- 0
           b <- 1
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- (a+b)/2
           VarX <- (b-a)^2/12
           q <- qunif(Ftheta,min=a,max=b)
         },
         Normal = {
           a <- dist.par[1]
           b <- dist.par[2]
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- a
           VarX <- b^2
           q <- qnorm(Ftheta, mean = a, sd = b)
         },
         Normal2 = {
           a <- dist.par[1]
           b <- dist.par[2]
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- a^2 + b^2
           VarX <- 4 * a^2 * b^2 + 2 * b^4
           q <- qchisq(Ftheta, 1)
         },
         DoubleExp = {
           a <- dist.par[1]
           b <- dist.par[2]
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- a
           VarX <- 2 * b^2
           #theta <- (qunif(Ftheta,1) - EX)/sqrt(VarX) * (sigma) + mu

         },
         DoubleExp2 = {
           a <- dist.par[1]
           b <- dist.par[2]
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- 2 * b^2 + a^2
           EY3 <- 6 * b^3 + 6 * a * b^2 + 5 * a^3 # Y is Laplace
           EY4 <- 24 * b^4 + 4 * a * EY3 - 6 * a^2 * (2 * b^2 + a^2) + 5 * a^4 # Y is Laplace
           VarX <- EY4 - EX^2

           #theta <- (qunif(Ftheta,1) - EX)/sqrt(VarX) * (sigma) + mu

         },
         LogNormal = {
           a <- dist.par[1] # logmean
           b <- dist.par[2] # logsigma
           if(!is.null(dist.par)){
             a <- dist.par[1]
             b <- dist.par[2]
           }
           EX <- exp(a + b^2 / 2)
           VarX <- exp(2 * (a + b^2)) - exp(2 * a + b^2)
           q <- qlnorm(Ftheta)
         },
         Gamma = {
           k <- dist.par[1] # beta in Casella
           o <- dist.par[2] # alpha in Casella
           if(!is.null(dist.par)){
             k <- dist.par[1]
             o <- dist.par[2]
           }
           EX <- k * o
           VarX <- o * k^2
           q <- qgamma(Ftheta, shape = o, scale = k)
         },
         Weibull = {
           k <- dist.par[1]
           l <- dist.par[2]
           if(!is.null(dist.par)){
             k <- dist.par[1]
             l <- dist.par[2]
           }
           EX <- l * gamma(1 + 1 / k)
           VarX <- l^2 * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2)
           q <- qweibull(Ftheta,shape = k, scale = l)
         },
         t = {
           v <- dist.par[1]
           if(!is.null(dist.par)){
             v <- dist.par[1]
           }
           EX <- 0
           VarX <- v/(v-2)

           q <- qt(Ftheta, v)
         },
         { # Normal (default)
           a <- dist.par[1]
           b <- dist.par[2]
           EX <- a
           VarX <- b^2
           q <- qnorm(Ftheta, mean = a, sd = b)
         }
  )

  theta <- (q - EX)/sqrt(VarX) * (sigma) + mu
  return(theta)
}



