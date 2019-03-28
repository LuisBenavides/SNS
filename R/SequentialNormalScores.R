#' @title Sequential normal scores
#' @description Transform a vector X into SNS using initial observations Y if available
#' SNS follow the order of X.
#' If ties, average ranks are used.
#' If Y = NULL, normal scores are set relative to X.
#'
#' @family aggregate functions
#' @seealso \code{\link{prod}} for products, \code{\link{cumsum}} for cumulative
#'   sums, and \code{\link{colSums}}/\code{\link{rowSums}} marginal sums over
#'   high-dimensional arrays.
#' @inheritParams NS
#' @param X.id is the id of the vector X
#' @export
#' @examples
#' # EXAMPLE CONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#' # [1] -0.52440051 -0.31863936  0.08964235
#'
#' # EXAMPLE CONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#' # [1] -0.52440051 -0.31863936  0.08964235
#'
#' # EXAMPLE UNCONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- NULL
#' Ftheta <- NULL
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#' # [1] -0.6045853 -0.3186394  0.0000000
#'
#' # EXAMPLE CONDITIONAL WITHOUT REFERENCE SAMPLE
#' Y <- NULL # c(10,20,30,40,50,60,70,80,90,100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#' # [1] -0.6744898 -0.3186394  0.6744898
#'
#' # EXAMPLE UNCONDITIONAL WITHOUT REFERENCE SAMPLE
#' Y <- NULL
#' X <- c(30, 35, 45)
#' theta <- NULL
#' Ftheta <- NULL
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#' # [1] 0.0000000 0.6744898 0.9674216
SNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL, alignment = "unadjusted", constant = NULL) {
  o.id <- unique(X.id) # originalid
  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  } else if (length(X) != length(X.id)) {
    print("ERROR, observations (X) have different length of the observations id (X.id)")
    return()
  }

  # detect the changes in the observation id vector
  changes.in.X.id <- c(1, as.numeric(X.id[1:(length(X.id) - 1)] != X.id[2:(length(X.id))]))
  # change the observation id
  X.id <- cumsum(changes.in.X.id)
  # get the different groups of the id
  groups <- unique(X.id)
  z <- rep(NA, length(groups)) # preallocate memory to initialize the SNS (one for group)
  i <- 1 # initialize the group index of the observation id vector
  if (is.null(Y)) { # if there is no reference sample
    Xe <- X[which(X.id == 1)] # get the first group
    ns <- NS(X = Xe, Y = NULL, theta = theta, Ftheta = Ftheta, alignment = alignment, constant = constant) # calculate the normal score
    z[i] <- sum(ns) / length(ns) # it is a vector with a subgroup size so it is needed to
    # make a mean
    Y <- Xe # the reference sample is the observations of the first group
    i <- i + 1
  }
  if (length(groups) > 1) { # if there are more than one group
    Ye <- Y # initialize the previous information to evaluate
    while (i <= length(groups)) { # repeat until the total groups are analized
      Xe <- X[which(X.id == groups[i])] # get the observations to evalute from the positions
      ad <- dataAlignment(Xe, Ye, alignment = alignment)
      Xe <- ad$X
      Ye <- ad$Y
      ns <- NS(X = Xe, Y = Ye, theta = theta, Ftheta = Ftheta, alignment = alignment, constant = constant) # calculate the normal score
      z[i] <- sum(ns) / length(ns) # it is a vector with a subgroup size so it is needed to
      # make a mean
      n <- length(Xe) # number of observations in the subgroup
      if (z[i] < 3 / sqrt(n) && z[i] > -3 / sqrt(n)) { # if the subgroup is in control
        Ye <- c(Ye, Xe) # add to reference sample the just evaluated observations
      }
      i <- i + 1 # continue with the next group
    }
  }

  return(z) # return the sequential normal score
}
