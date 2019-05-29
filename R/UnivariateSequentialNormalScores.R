#' @title Sequential Normal Scores
#' @description Transform a vector \code{X} into SNS using initial observations \code{Y} if available
#' SNS follow the order of \code{X}.
#' @section Comments:
#' If ties, average ranks are used.
#' @seealso \code{\link{NS}} for normal scores
#' @inheritParams NS
#' @inheritParams getRL
#' @param X.id vector. The id of the vector \code{X}.
#' @param isFixed logical. If \code{TRUE} the reference sample does not update, otherwise the reference sample is updated when the batch is in control.
#' @param snsRaw logical. If \code{TRUE} return also the sns for each observation in vector \code{X}.
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
SNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL, scoring = "Z",
                alignment = "unadjusted", constant = NULL, absolute = FALSE,
                chart="Shewhart", chart.par=c(3),
                snsRaw = FALSE, isFixed = FALSE) {

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
  if(snsRaw){
    Zraw = rep(NA, length(X))
  }
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

  ucl = 0
  if (scoring == "Z-SQ"){
    inf = 1000000 #infinite value to better approximation
    alpha = 0.005 #confindent interval
    vec <- rchisq(inf, length(X[which(Xb.id == groups[i])])) #chi-sq random generator numbers according to the "infinite value"
    ucl <- quantile(vec , 1-alpha) #control limit
  }
  while (i <= length(groups)) { # repeat until the total groups are analized
    Xb = X[which(Xb.id == groups[i])] # get the observations to evalute from the positions
    ad = SNS::dataAlignment(Xb, Yb, alignment = alignment)
    Xb = ad$X
    Yb = ad$Y
    ns = SNS::NS(X = Xb, Y = Yb, theta = theta, Ftheta = Ftheta, scoring = scoring, alignment = alignment, constant = constant) # calculate the normal score
    ns = ns$Z

    n = length(Xb)
    if(snsRaw){
      Zraw[(1+n*(i-1)):(n+n*(i-1))] = ns
    }

    switch (scoring,
      "Z" = {# it is a vector with a subgroup size so it is needed to average them
        z[i] = mean(ns)
      },
      "Z-SQ" = {# it is a vector with a subgroup size so it is needed to sum them
        z[i] = sum(ns)
      }
    )
    Z = z[i]
    if (is.null(Yb) && i == 1) { # if there is no reference sample
      Yb = Xb
    }
    # check if the subgroup is in control according to each scheme
    # the reference sample is updated
    updateSample <- FALSE
    switch(chart,
           Shewhart = {
             if (scoring == "Z"){
               ucl = k / sqrt(n)
             }
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
             Cminus[i] <- cminus

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
    if (scoring == "Z-SQ"){
      LCL[i] = 0
    }


    if (updateSample && !isFixed){# if the subgroup is in control (updateSample change to TRUE)
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
    LCL = LCL,
    scoring = scoring
  )
  if(snsRaw){
    output$Zraw = Zraw
  }
  switch(chart,
         CUSUM = {
           output$Cplus = Cplus
           output$Cminus = Cminus
         },
         EWMA = {
           output$E = E
         }
  )
  class(output)="SNS" # Class definition
  return(output) # return the sequential normal score
}

#' @import graphics
#' @export
plot.SNS <- function(x,...){
  par(mar = c(6,6,4,2))

  Z = x$Z
  n = x$n
  o.id = unique(x$X.id) # original id
  chart = coef(x)$chart
  chart.par = coef(x)$chart.par
  UCL = x$UCL
  LCL = x$LCL
  Cplus = x$Cplus
  Cminus = x$Cminus
  E = x$E
  scoring = x$scoring
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
           ylab = "Z"
          if (scoring == "Z-SQ"){
            ymin = 0
            ylab = expression(Z^2)
          }
           plot(o.id, Z, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab=ylab,cex.lab=2.5, cex.axis=1.5, cex=2)
           lines(o.id, Z, lt=2, lwd=3)
         },
         CUSUM = {
           type = chart.par[3]
           switch(type,
                  "1" = {
                    plot(o.id, Cplus, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab=expression(C^"+"),cex.lab=2.5, cex.axis=1.5, cex=2)
                    lines(o.id, Cplus, lt=2, lwd=3)
                  },
                  "2" = {
                    plot(o.id, Cminus, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab=expression(C^"-"),cex.lab=2.5, cex.axis=1.5, cex=2)
                    lines(o.id, Cminus, lt=2, lwd=3)
                  },
                  "3" = {
                    plot(o.id, Cplus, pch=19, ylim=c(ymin, ymax), xlab = "Batch",ylab=expression(C^"+"/C^"-"),cex.lab=2.5, cex.axis=1.5, cex=2)
                    points(o.id, Cminus, pch=15, cex=2)
                    lines(o.id, Cplus, lt=2, lwd=3)
                    lines(o.id, Cminus, lt=2, lwd=3)
                    legend("topleft", c(expression(C^"+"), expression(C^"-")), pch=c(19, 15))
                  }
           )
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
  if(sum(LCL) != 0)  lines(o.id, LCL,lt=4, lwd=3)
}