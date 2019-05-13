#' @title Multivariate Run Length
#' @description Get the run length
#' @inheritParams mGetDist
#' @inheritParams MNS
#' @param replica scalar. It is used for the parallel version of the function (\code{parallel=TRUE}). Default \code{1}.
#' @param n scalar. Subroup size
#' @param m scalar. Reference sample size
#' @param chart character string. Selected type of chart. One option available: T2
#' \describe{
#'   \item{T2 scheme: }{is \code{c(k)}, where \code{k} comes from \eqn{UCL = mu + k\sigma, LCL = mu - k\sigma.}}
#' }
#' @param chart.par vector. Control limit and other parameters of the selected chart.
#' @param calibrate logical. If \code{TRUE} the RL is limit to 10 times the target ARL.
#' @param arl0 scalar. Expected value of the RL. Default \code{370}.
#' @export
#' @import stats
#' @examples
#' nv <- 2 #number of variables
#' inf <- 1000000 #infinite value to better approximation
#' alpha <- 0.005 #confindent interval
#' vec <- rchisq(inf, nv) #chi-sq random generator numbers according to the "infinite value"
#' h <- quantile(vec , 1-alpha) #control limit
#' mGetRL(n=5, m=10, nv=nv, mu=c(0,0), dists = c("Normal", "Normal"),
#' dists.par = matrix(c(0,1,1,0,1,1), ncol=2),
#' chart = "T2",chart.par=c(h), correlation=0)
mGetRL <- function(replica = 1, n, m, nv, theta = NULL, Ftheta = NULL,
                  dists, mu, sigma=NULL, dists.par = NULL, correlation=0,
                  chart="T2", chart.par = c(10),
                  alignment = "unadjusted", constant = NULL, absolute=FALSE,
                  calibrate=FALSE, arl0=370) {
  # initilize the reference sample
  Y <- NULL
  Z = NULL #preallocate memory for sns

  if (m > 0) { # if there are reference sample
    # generate the reference sample
    Y <- mGetDist(n = m, nv = nv, mu = mu[1], sigma = sigma, correlation=correlation, dists = dists, dists.par = dists.par)
    ns <- MNS(X = Y, Y = NULL, theta = theta, Ftheta = Ftheta, alignment = alignment, constant = constant)
    Z <- ns$Z
  }

  RL <- 0
  in.Control <- TRUE
  switch(chart,
     T2 = {
       ucl <- chart.par[1] #control limit
     }
  )

  while (in.Control) {
    # add one iteration to run length
    RL <- RL + 1

    # generate the subgroup to monitor
    X <- mGetDist(n = n, nv = nv, mu = mu[2],sigma=sigma, dists = dists, dists.par = dists.par, correlation=correlation)

    # get the normal scores
    ns <- MNS(X = X, Y = Y, theta = theta, Ftheta = Ftheta, alignment = alignment, constant = constant)
    Zb <- ns$Z


    if (is.null(Y)) { # if is the first batch
      T2 = 0 # it does not give any information and is considered the reference sample
    }else{

      muZ <- apply(Zb, 2, mean) #obtain the mean for each variable in the batch
      # check if the subgroup is in control according to each scheme
      # the reference sample is updated

      T2 <- n*(muZ%*%chol2inv(chol(cor(Z, method = "spearman")))%*%muZ) #get the T2 statistic
    }


    # if the subgroup is out of the limits
    # an alarm is detected
    switch(chart,
       T2 = {
         # if the subgroup is out of the limits an alarm is detected
         if (T2 > ucl) in.Control <- FALSE
       }
    )
    #cat("RL=",RL, " T2=",T2," ucl",ucl, "\n")
    if (calibrate) if (RL >= arl0 * 50) in.Control <- FALSE
    if (RL >= arl0 * 1000) in.Control <- FALSE

    Y <- rbind(Y, X) # update the reference sample
    Z <- rbind(Z, Zb) #update the sns (for correlation matrix)
  }
  return(RL)
}


#' @title Multivariate Average Run Length (ARL)
#' @description Get the ARL \code{\link{getRL}}
#' @inheritParams mGetRL
#' @param print.RL logical. If \code{TRUE} return the vectors of RL for each iteration.
#' @param replicates scalar. Number of replicates to get the ARL
#' @param progress logical. If \code{TRUE} it shows the progress in the console.
#' @param isParallel logical. If \code{TRUE} the code runs in parallel according to the
#' number of cores in the computer,otherwise the code runs sequentially. Default \code{TRUE}.
#' @export
#' @import parallel
#' @import stats
#' @examples
#' mGetARL(replicates=50,n=5,m=100,nv=2,mu=c(0,0),
#' dists = c("Normal", "Normal"), dists.par = matrix(c(0,1,1,0,1,1), ncol=2),
#' isParallel=FALSE)
mGetARL <- function(n, m, nv, theta = NULL, Ftheta = NULL,
                   dists, dists.par = NULL, mu, sigma=NULL,
                   chart = "T2", chart.par = c(10), correlation = 0,
                   replicates = 10000, isParallel = TRUE,
                   print.RL = FALSE, progress = FALSE,
                   calibrate = FALSE, arl0 = 370,
                   alignment = "unadjusted", constant = NULL, absolute=FALSE) {
  RLs <- NULL
  if (isParallel) {
    cluster <- makeCluster(detectCores() - 1)
    clusterExport(cluster, "MNS")
    clusterExport(cluster, "mGetDist")
    clusterExport(cluster, "mGetRL")
    RLs <- parSapply(cluster, 1:replicates, mGetRL, n = n, m = m, nv = nv, theta = theta, Ftheta = Ftheta, dists = dists, mu = mu, dists.par = dists.par, chart = chart, chart.par=chart.par,correlation=correlation, alignment=alignment, constant=constant,absolute=absolute)
    stopCluster(cluster)
  } else {
    t0 <- Sys.time()
    for (r in 1:replicates) {
      RL <- mGetRL(1, n = n, m = m, nv = nv, theta = theta, Ftheta = Ftheta, dists = dists, mu = mu, dists.par = dists.par, chart = chart, alignment=alignment, constant=constant,absolute=absolute)

      RLs <- c(RLs, RL)

      # print out progress
      if (progress) { # if is TRUE
        if (r %% 10 == 0) { # every 10 replicates
          t1 <- Sys.time()
          remaining.iterations <- replicates - r
          remaining.time <- remaining.iterations * difftime(t1, t0, units = "min") / r
          cat("ARL", round(mean(RLs), digits = 1), "-- SDRL", round(sd(RLs), digits = 1), "--> Time remaining", remaining.time, "in minutes to complete", remaining.iterations, "iterations", "\n", sep = " ")
        }
      }
    }
  }

  output <- list(
    ARL = mean(RLs),
    SDRL = sd(RLs),
    MRL = median(RLs),
    QRL = quantile(x = RLs, probs = c(0.05, 0.1, 0.2, 0.25, 0.5, 0.75, 0.8, 0.9, 0.95), names = TRUE, type = 3)
  )
  if (print.RL) output$RL <- RLs

  if (progress) cat("Final ARL", round(mean(RLs), digits = 1), "-- SDRL", round(sd(RLs), digits = 1), "\n", "See output variable for more.\n\n", sep = " ")

  return(output)
}
