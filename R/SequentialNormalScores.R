#' @title Sequential Normal Scores
#' @description Transform a vector \code{X} into SNS using initial observations \code{Y} if available
#' SNS follow the order of \code{X}.
#' @section Comments:
#' If ties, average ranks are used.
#' @seealso \code{\link{NS}} for normal scores
#' @inheritParams NS
#' @param X.id vector. The id of the vector \code{X}.
#' @export
#' @examples
#' # EXAMPLE CONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#'
#' # EXAMPLE CONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#'
#' # EXAMPLE UNCONDITIONAL WITH REFERENCE SAMPLE
#' Y <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
#' X <- c(30, 35, 45)
#' theta <- NULL
#' Ftheta <- NULL
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#'
#' # EXAMPLE CONDITIONAL WITHOUT REFERENCE SAMPLE
#' Y <- NULL # c(10,20,30,40,50,60,70,80,90,100)
#' X <- c(30, 35, 45)
#' theta <- 40
#' Ftheta <- 0.5
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#'
#' # EXAMPLE UNCONDITIONAL WITHOUT REFERENCE SAMPLE
#' Y <- NULL
#' X <- c(30, 35, 45)
#' theta <- NULL
#' Ftheta <- NULL
#' sample.id <- c("a", "b", "c")
#' SNS(X = X, X.id = sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
SNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL,
                alignment = "unadjusted", constant = NULL, absolute = FALSE,
                chart="Shewhart", chart.par=c(3)) {

  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  } else if (length(X) != length(X.id)) {
    print("ERROR, observations (X) have different length of the observations id (X.id)")
    return()
  }

  # detect the changes in the observation id vector
  changes.in.X.id = c(1, as.numeric(X.id[1:(length(X.id) - 1)] != X.id[2:(length(X.id))]))
  #change the observation id
  Xb.id = cumsum(changes.in.X.id)
  #get the different groups of the id
  groups = unique(Xb.id)
  z = rep(NA, length(groups)) # preallocate memory to initialize the SNS (one for group)
  i = 1 # initialize the group index of the observation id vector
  Yb = Y
  if(!is.null(Yb)){
    Yb = Yb[!is.na(Yb)] # initialize reference sample (remove na values)
  }

  UCL = rep(NA, length(groups))
  LCL = rep(NA, length(groups))
  switch(chart,
         Shewhart = {
           k <- chart.par[1]
         },
         CUSUM = {
           k <- chart.par[1]
           h <- chart.par[2]
           type <- chart.par[3]
           Cplus = rep(NA, length(groups))
           Cminus = rep(NA, length(groups))
           cplus <- 0
           cminus <- 0
         },
         EWMA = {
           lambda <- chart.par[1]
           L <- chart.par[2]
           E = rep(NA, length(groups))
           e <- 0
         }
  )

  while (i <= length(groups)) { # repeat until the total groups are analized
    Xb = X[which(Xb.id == groups[i])] # get the observations to evalute from the positions
    ad = dataAlignment(Xb, Yb, alignment = alignment)
    Xb = ad$X
    Yb = ad$Y
    ns = NS(X = Xb, Y = Yb, theta = theta, Ftheta = Ftheta, alignment = alignment, constant = constant) # calculate the normal score
    ns = ns$Z
    n = length(Xb)
    z[i] = sum(ns) / n # it is a vector with a subgroup size so it is needed to
    Z = z[i]
    # make a mean
    if (is.null(Yb) && i == 1) { # if there is no reference sample
      Yb = Xb
    }
    # check if the subgroup is in control according to each scheme
    # the reference sample is updated
    updateSample <- FALSE
    switch(chart,
           Shewhart = {
             ucl = k / sqrt(n)
             if (abs(Z) < ucl) updateSample <- TRUE
           },
           CUSUM = {
             switch(type,
                    "1" = {
                      cplus <- max(c(0, cplus + Z * sqrt(n) - k))
                    },
                    "2" = {
                      cminus <- min(c(0, cminus + Z * sqrt(n) + k))
                    },
                    "3" = {
                      cplus <- max(c(0, cplus + Z * sqrt(n) - k))
                      cminus <- min(c(0, cminus + Z * sqrt(n) + k))
                    }
             )

             Cplus[i] <- cplus
             Cminus[i] = cminus

             ucl = h

             if (cplus < ucl || cminus > -ucl) updateSample <- TRUE
           },
           EWMA = {
             e <- lambda * Z + (1 - lambda) * e

             E[i] = e

             ucl <- L / sqrt(n) * sqrt(lambda / (2 - lambda) * (1 - (1 - lambda)^(2 * i)))
             if (abs(e) < ucl) updateSample <- TRUE
           }
    )

    UCL[i] = ucl
    LCL[i] = -ucl

    if (updateSample){# if the subgroup is in control (updateSample change to TRUE)
      Yb = c(Yb, Xb) # add to reference sample the new observations
    }
    i = i + 1 # continue with the next group
  }
  output = list(
    coefficients = list(
      n=n,
      chart = chart,
      chart.par = chart.par
    ),
    Z = z,
    X.id = X.id,
    UCL = UCL,
    LCL = LCL
  )
  switch(chart,
         CUSUM = {
           output$Cplus = Cplus
           output$Cminus = Cminus
         },
         EWMA = {
           output$E = E
         }
  )

  class(output)="sns" # Class definition

  return(output) # return the sequential normal score
}


#' @title Plot Sequential Normal Scores
#' @description plot the Sequential Normal Scores by using only \code{plot}
#' @import graphics
plot.sns=function(object,...){
  par(mar = c(6,6,4,2))

  Z = object$Z
  n = object$n
  o.id = unique(object$X.id) # original id
  chart = coef(object)$chart
  chart.par = coef(object)$chart.par
  UCL = object$UCL
  LCL = object$LCL
  Cplus = object$Cplus
  Cminus = object$Cminus
  E = object$E
  switch(chart,
         Shewhart = {
           difMaxZ = abs(max(Z) - max(UCL))
           difMinZ = abs(min(Z) - min(LCL))
         },
         CUSUM = {
           difMaxZ = abs(max(c(Cplus, Cminus)) - max(UCL))
           difMinZ = abs(min(c(Cplus, Cminus)) - min(LCL))
         },
         EWMA = {
           difMaxZ = abs(max(E) - max(UCL))
           difMinZ = abs(min(E) - min(LCL))
         }
  )

  ymax = max(UCL)
  if(difMaxZ > difMinZ){
    ymax = ymax + difMaxZ
  }else{
    ymax = ymax + difMinZ
  }

  ymin = min(LCL)
  if(difMaxZ < difMinZ){
    ymin = ymin - difMaxZ
  }else{
    ymin = ymin - difMinZ
  }

  switch(chart,
         Shewhart = {
           plot(o.id, Z, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab="Z",cex.lab=2.5, cex.axis=1.5, cex=2)
           lines(o.id, Z, lt=2, lwd=3)
         },
         CUSUM = {
           plot(o.id, Cplus, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab="Cplus/Cminus",cex.lab=2.5, cex.axis=1.5, cex=2)
           points(o.id, Cminus, pch=15, cex=2)
           lines(o.id, Cplus, lt=2, lwd=3)
           lines(o.id, Cminus, lt=2, lwd=3)
         },
         EWMA = {
           plot(o.id, E, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab="E",cex.lab=2.5, cex.axis=1.5, cex=2)
           lines(o.id, E, lt=2, lwd=3)
         }
  )


  change = 1
  if(o.id[1] > o.id[length(o.id)]){
    change = -1
  }
  lines(o.id, UCL,lt=4, lwd=3)
  lines(o.id, LCL,lt=4, lwd=3)
  #lines(c(o.id[1]-10*change,o.id,o.id[length(o.id)]+10*change), rep(UCL,length(o.id)+2),lt=4, lwd=3)
  #lines(c(o.id[1]-10*change,o.id,o.id[length(o.id)]+10*change), rep(LCL,length(o.id)+2),lt=4, lwd=3)
}
