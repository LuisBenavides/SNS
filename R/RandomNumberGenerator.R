#' @title Random Observations Generetor
#' @description Random observations generator selected from several distributions with user defined mean and variance.
#' @param n scalar. Number of observations to be generated.
#' @param dist character string. Select from:
#' \itemize{
#'   \item{"Normal": Normal distribution (default).}
#'   \item{"DoubleExp": Double exponential distribution (also known as Laplace distribution).}
#'   \item{"DoubleExp2": Double exponential squared distribution from a \code{DoubleExp(0,1)}.}
#'   \item{"LogNormal": Lognormal distribution.}
#'   \item{"Gamma": Gamma distribution.}
#'   \item{"Weibull": Weibull distribution.}
#' }
#' @param mu scalar. Expected value of the desired distribution.
#' @param sigma scalar. Standard deviation of the desired distribution.
#' @param par.location scalar. Location parameter of the desired distribution. Default 0**.
#' @param par.scale scalar. Scale parameter of the desired distribution. Default 1**.
#' @param par.shape scalar. Shape parameter of the desired distribution, Default 1.
#' @section **Note:
#' \itemize{
#'   \item{For "Lognormal", \code{par.location} and \code{par.scale} correspond to the location and scale parameters of the normal
#'     distribution that generales the lognormal. Hence, in this case they are the logmean and
#'     the logsigma parameters}
#'   \item{For "Normal2" and "DoubleExp2", \code{par.location} and \code{par.scale} correspond
#'     correspond to the location and scale parameters of the normal and double exponential
#'     that are used to generates their squared forms.}
#' }
#' @return A vector \code{x} with \code{n} observations generated following the selected distribution with its parameters.
#' @export
#' @examples
#' getDist(1, "Normal", 0, 1)
getDist <- function(n, dist, mu, sigma,
                    par.location = 0, par.scale = 1, par.shape = 1) {
  switch(dist,
    Normal = {
      a <- par.location
      b <- par.scale
      EX <- a
      VarX <- b^2
      z <- (rnorm(n, mean = a, sd = b) - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    Normal2 = {
      a <- par.location
      b <- par.scale
      EX <- a^2 + b^2
      VarX <- 4 * a^2 * b^2 + 2 * b^4
      z <- ((rnorm(n, mean = a, sd = b))^2 - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    DoubleExp = {
      a <- par.location
      b <- par.scale
      EX <- a
      VarX <- 2 * b^2
      # xtemp = a - b * sign(U - 0.5) * log(1 - 2 * abs(U - 0.5)) #This one appeared in Wikipedia
      xtemp <- log(runif(n) / runif(n)) / 2^(0.5) # this is the recommended method. Gives standard DE variates.

      z <- (xtemp - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    DoubleExp2 = {
      a <- par.location
      b <- par.scale
      EX <- 2 * b^2 + a^2
      EY3 <- 6 * b^3 + 6 * a * b^2 + 5 * a^3 # Y is Laplace
      EY4 <- 24 * b^4 + 4 * a * EY3 - 6 * a^2 * (2 * b^2 + a^2) + 5 * a^4 # Y is Laplace
      VarX <- EY4 - EX^2
      xtemp <- (log(runif(n) / runif(n)))^2
      z <- (xtemp - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    LogNormal = {
      a <- par.location # logmean
      b <- par.scale # logsigma
      EX <- exp(a + b^2 / 2)
      VarX <- exp(2 * (a + b^2)) - exp(2 * a + b^2)
      xtemp <- rlnorm(n, meanlog = a, sdlog = b)
      # xtemp = exp(a + b*rnorm(n))
      z <- (xtemp - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    Gamma = {
      k <- par.scale # beta in Casella
      o <- par.shape # alpha in Casella
      EX <- k * o
      VarX <- o * k^2
      xtemp <- rgamma(n, shape = o, scale = k)
      z <- (xtemp - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    Weibull = {
      k <- par.shape
      l <- par.scale
      EX <- l * gamma(1 + 1 / k)
      VarX <- l^2 * (gamma(1 + 2 / k) - (gamma(1 + 1 / k))^2)
      xtemp <- rweibull(n, shape = k, scale = l)
      z <- (xtemp - EX) / VarX^(0.5)
      x <- mu + sigma * z
    },
    foo = {

    }, { # Normal (default)
      a <- par.location
      b <- par.scale
      EX <- a
      VarX <- b^2
      z <- (rnorm(n, mean = a, sd = b) - EX) / VarX^(0.5)
      x <- mu + sigma * z
    }
  )
  return(x)
}
