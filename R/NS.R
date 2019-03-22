#' Normal scores X from Y
#'
#' Get conditional and unconditional normal score (NS) of observations (X)
#' relative to previous observations (Y).
#'
#' If Y = NULL (no previous observation available), NS is relative to X.
#' GENERAL COMMENTS
#' If ties, average ranks are used.
#' If Y = NULL, normal scores are set relative to X.
#' UNCONDITIONAL COMMENTS
#' Instead of Van Der Waerden Normal Scores where p = r/(n+1), p = (r-0.5)/n,
#' where r stands for rank and p for the input evaluated in the
#' inverse of a Standard Normal Distribution.
#'
#' @param X is a numerical vector.
#' @param Y is a numerical vector. If Y is not defined, by default Y = NULL.
#' @param theta is a constant.
#' @param Ftheta is a constant between (0,1).
#'
#' @return Z (Normal Scores)
#' @export
#' @examples
#' Y = c(10,20,30,40,50,60,70,80,90,100)
#' X = c(30, 35, 45)
#' theta = 40
#' Ftheta = 0.5
#' #EXAMPLE CONDITIONAL
#' NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
#' #[1] -0.52440051 -0.38532047  0.08964235
#' #EXAMPLE UNCONDITIONAL
#' theta = NULL
#' Ftheta = NULL
#' NS(X = X, Y = Y, theta = theta, Ftheta = Ftheta)
#' #[1] -0.6045853 -0.4727891 -0.2298841

NS <- function(X, Y = NULL, theta = NULL, Ftheta = NULL){
  #Check for errors
  if(is.null(theta) != is.null(Ftheta)){ #in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  }
  n = length(X)       #get the number of observations
  if(is.null(theta) | is.null(Ftheta)){ #if descriptive data is not vailable
                                              #such as a quantil (theta) or
    if(is.null(Y)){   #if previous data is not available
      r = rank(X)     #rank the observations with general wanking function
    }else{            #if previous data is available
      r = rep(NA, n)  #preallocate memory to initialize the ranks. One for each observation.
      for(i in 1:n){  #for each observation, by index.
        #obtain the rank by comparing each obsarvation
        #depending on if is greater or equals to previous data
        r[i] = sum(Y < X[i]) + (sum(Y == X[i]) + 2)/2
      }
      n = length(Y) + 1 #uptade number of observations and add one unit
    }
    p = (r-0.5)/n       #obtain the probability of the ranks
  }else{
    if(is.null(Y)){ #if previous data is not available
      Y = X         #previous data is the observed data
    }
    #Nminus = sum(Y <= theta) #numbers of <= theta used in individual ranking
    #Nplus = sum(Y > theta) #number > theta used in individual ranking.
    nX = length(X)  #obtain the number of normal scores needed
    r = rep(NA, n)  #preallocate memory to initialize the ranks. One for each observation.
    p = rep(NA, n)  #preallocate memory to initialize the probability. One for each observation.
    for(i in 1:n){  #for each observation, by index.
      r[i] = (sum(Y < X[i] & Y <= theta) + (sum(Y == X[i] & Y <= theta) + 2)/2)* (X[i] <= theta) + (sum(Y < X[i] & Y > theta) + (sum(Y == X[i] & Y > theta) + 2)/2) * (X[i] > theta)
      nTheta = (X[i] <= theta) * sum(Y <= theta) + (X[i] > theta) * sum(Y > theta) + 1
      p[i] = Ftheta * (X[i] > theta) + ( (1 - Ftheta) * (X[i] > theta) + (X[i] <= theta) * Ftheta ) * (r[i] - 0.5) / nTheta
    }
  }
  z = qnorm(p)      #evaluated the inverse of a Standard Normal Distribution of the probability
                    #to obtain the Normal Scores (NS).
  return(z)
}


