#' @title Multivariate Normal Scores
#' @description Get conditional or unconditional multivariate normal score (NS) of observations (\code{X})
#' relative to previous observations (\code{Y}).
#' @inheritParams dataAlignment
#' @param X matrix or data.frame. New observations to obtain the normal scores.
#' @param Y matrix or data.frame. If \code{Y} is not defined (no previous observation available, \code{NULL}), NS is relative to \code{X}. Default \code{NULL}.
#' @param theta vector. Value corresponding with the \code{Ftheta} quantile.
#' @param Ftheta vector. Quantile of the data distribution. The values that take are between (0,1).
#' @param scoring character string. If "Z" (normal scores) (default). If "Z-SQ" (normal scores squared).
#' @return Multiple output. Select by \code{output$}
#' \itemize{
#'   \item \code{R}: matrix. Multivariate Ranks for the \code{X} observations. If ties occurs, average ranks are used.
#'   \item \code{P}: matrix. Multivariate Probability of the ranks for the \code{X} observations. Instead of Van Der Waerden normal scores where \eqn{P = R/(n+1)}, \eqn{P = (R-0.5)/n},
#' where \eqn{R} stands for rank and \eqn{P} for the input evaluated in the inverse of a Standard Normal Distribution.
#'   \item \code{Z}: matrix. Multivariate Normal scores for the \code{X} observations. \eqn{Z} if \code{scoring} is "Z" and \eqn{Z^2} if \code{scoring} is "Z-SQ".
#' }
#' @export
#' @examples
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' Y = matrix(Y, ncol=2)
#' X <- c(30, 35, 45, 30, 35, 45)
#' X = matrix(X, ncol=2)
#' theta <- c(40, 40)
#' Ftheta <- c(0.5, 0.5)
#' # EXAMPLE CONDITIONAL
#' MNS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)

MNS <- function(X, Y = NULL, theta = NULL, Ftheta = NULL, scoring = "Z",
               alignment = "unadjusted", constant = NULL, absolute = FALSE) {
  # Check for errors
  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  }

  ad <- SNS::dataAlignment(X=X, Y=Y, alignment = alignment, constant = constant) # Alignment of the data
  X <- ad$X
  Y <- ad$Y

  nv = ncol(X) #get the number of variables
  n <- nrow(X) # get the number of observations
  if (is.null(theta) | is.null(Ftheta)) { # if descriptive data is not available
    # such as a quantil (theta)
    if (is.null(Y)) { # if previous data is not available
      R <- apply(X, 2, rank) # rank the observations with general wanking function
    } else { # if previous data is available
      R = matrix(0, n, nv) #preallocate memory to initialize the ranks.
      for(i in 1:n){#for each observation in the batch
        for (j in 1:nv){#for each variable
          #for each observation, by index.
          #obtain the rank by comparing each observation
          #depending on if is greater or equals to previous data
          R[i,j] = sum(Y[,j] < X[i,j]) + (sum(Y[,j] == X[i,j]) + 2)/2
        }
      }
      n <- nrow(Y) + 1 # update number of observations and add one unit
    }
    P <- (R - 0.5) / n # obtain the probability of the ranks
  } else {
    if (is.null(Y)) { # if previous data is not available
      Y <- X # previous data is the observed data
    }
    # Nminus = sum(Y <= theta) #numbers of <= theta used in individual ranking
    # Nplus = sum(Y > theta) #number > theta used in individual ranking.
    R = matrix(0, n, nv) # preallocate memory to initialize the ranks. One for each observation.
    P = matrix(0, n, nv) # preallocate memory to initialize the probability. One for each observation.
    for (i in 1:n) { # for each observation, by index.
      for (j in 1:nv) { # for each variable, by index.
        R[i,j] <- (sum(Y[,j] < X[i,j] & Y[,j] <= theta[j]) + (sum(Y[,j] == X[i,j] & Y[,j] <= theta[j]) + 2) / 2) * (X[i,j] <= theta[j]) + (sum(Y[,j] < X[i,j] & Y[,j] > theta[j]) + (sum(Y[,j] == X[i,j] & Y[,j] > theta[j]) + 2) / 2) * (X[i,j] > theta[j])
        nTheta <- (X[i,j] <= theta[j]) * sum(Y[,j] <= theta[j]) + (X[i,j] > theta[j]) * sum(Y[,j] > theta[j]) + 1
        P[i,j] <- Ftheta[j] * (X[i,j] > theta[j]) + ((1 - Ftheta[j]) * (X[i,j] > theta[j]) + (X[i,j] <= theta[j]) * Ftheta[j]) * (R[i,j] - 0.5) / nTheta
      }
    }
  }

  Z <- qnorm(P) # evaluated the inverse of a Standard Normal Distribution of the probability
  # to obtain the Multivariate Normal Scores (MNS).

  switch(scoring,
         "Z-SQ" = {
           Z <- Z^2
         },
         "Z" = {
           Z <- Z
         },
         {#default
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
