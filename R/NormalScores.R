#' @title Normal Scores
#' @description Get conditional or unconditional normal score (NS) of observations (\code{X})
#' relative to previous observations (\code{Y}).
#' @inheritParams dataAlignment
#' @param X vector. New observations to obtain the N¡normal scores.
#' @param Y vector. If \code{Y} is not defined (no previous observation available, \code{NULL}), NS is relative to \code{X}. Default \code{NULL}.
#' @param theta scalar. Value corresponig with the \code{Ftheta} quantile.
#' @param Ftheta scalar. Quantile of the data distribution. The values that take are between (0,1).
#' @param scoring character string. If "Z" (normal scores) (default). If "Z-SQ" (normal scores squared).
#' @return Multiple output. Select by \code{output$}
#' \itemize{
#'   \item \code{R}: vector. Ranks for the \code{X} observations. If ties occurs, average ranks are used.
#'   \item \code{P}: vector. Probability of the ranks for the \code{X} observations. Instead of Van Der Waerden normal scores where \eqn{P = R/(n+1)}, \eqn{P = (R-0.5)/n},
#' where \eqn{R} stands for rank and \eqn{P} for the input evaluated in the inverse of a Standard Normal Distribution.
#'   \item \code{Z}: vector. Normal scores for the \code{X} observations. \eqn{Z} if \code{scoring} is "Z" and \eqn{Z^2} if \code{scoring} is "Z-SQ".
#' }
#' @export
#' @examples
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' # EXAMPLE CONDITIONAL
#' NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
#'
#' # EXAMPLE UNCONDITIONAL
#' theta <- NULL
#' Ftheta <- NULL
#' NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
NS <- function(X, Y = NULL, theta = NULL, Ftheta = NULL, scoring = "Z",
               alignment = "unadjusted", constant = NULL, absolute = FALSE) {
  # Check for errors
  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  }

  ad <- SNS::dataAlignment(X, Y, alignment = alignment, constant = constant) # Alignment of the data
  X <- ad$X
  Y <- ad$Y

  n <- length(X) # get the number of observations
  if (is.null(theta) | is.null(Ftheta)) { # if descriptive data is not vailable
    # such as a quantil (theta) or
    if (is.null(Y)) { # if previous data is not available
      R <- rank(X) # rank the observations with general wanking function
    } else { # if previous data is available
      R <- rep(NA, n) # preallocate memory to initialize the ranks. One for each observation.
      for (i in 1:n) { # for each observation, by index.
        # obtain the rank by comparing each obsarvation
        # depending on if is greater or equals to previous data
        R[i] <- sum(Y < X[i]) + (sum(Y == X[i]) + 2) / 2
      }
      n <- length(Y) + 1 # uptade number of observations and add one unit
    }
    P <- (R - 0.5) / n # obtain the probability of the ranks
  } else {
    if (is.null(Y)) { # if previous data is not available
      Y <- X # previous data is the observed data
    }
    # Nminus = sum(Y <= theta) #numbers of <= theta used in individual ranking
    # Nplus = sum(Y > theta) #number > theta used in individual ranking.
    R <- rep(NA, n) # preallocate memory to initialize the ranks. One for each observation.
    P <- rep(NA, n) # preallocate memory to initialize the probability. One for each observation.
    for (i in 1:n) { # for each observation, by index.
      R[i] <- (sum(Y < X[i] & Y <= theta) + (sum(Y == X[i] & Y <= theta) + 2) / 2) * (X[i] <= theta) + (sum(Y < X[i] & Y > theta) + (sum(Y == X[i] & Y > theta) + 2) / 2) * (X[i] > theta)
      nTheta <- (X[i] <= theta) * sum(Y <= theta) + (X[i] > theta) * sum(Y > theta) + 1
      P[i] <- Ftheta * (X[i] > theta) + ((1 - Ftheta) * (X[i] > theta) + (X[i] <= theta) * Ftheta) * (R[i] - 0.5) / nTheta
    }
  }
  Z <- qnorm(P) # evaluated the inverse of a Standard Normal Distribution of the probability
  # to obtain the Normal Scores (NS).

  switch(scoring,
    "Z-SQ" = {
      Z <- Z^2
    },
    "Z" = {
      Z <- Z
    },
    {
      Z = Z
    }
  )

  output <- list(
    R = R,
    P = P,
    Z = Z
  )
  return(output)
}

