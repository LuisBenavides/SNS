#' @title Multivariate Sequential Normal Scores
#' @description Transform a matrix \code{X} into SNS using initial observations \code{Y} if available
#' SNS follow the order of \code{X}.
#' @section Comments:
#' If ties, average ranks are used.
#' @seealso \code{\link{MNS}} for multivariate normal scores
#' @inheritParams MNS
#' @inheritParams mgetRL
#' @param X.id vector. The id of each column (variable) of the matrix \code{X}.
#' @param isFixed logical. If \code{TRUE} the reference sample does not update, otherwise the reference sample is updated when the batch is in control.
#' @export
#' @examples
#' X = cbind(example91$X1, example91$X2)
#' X.id = example91$X1.id
#' msns = MSNS(X, X.id)
MSNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL, scoring = "Z",
                alignment = "unadjusted", constant = NULL, absolute = FALSE,
                chart="T2", chart.par = c(0.005), null.dist="Chi", isFixed = FALSE) {

  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  } else if (nrow(X) != length(X.id)) {
    print("ERROR, observations (X) have different length of the observations id (X.id)")
    return()
  }

  # detect the changes in the observation id vector
  changes.in.X.id = c(1, as.numeric(X.id[1:(length(X.id) - 1)] != X.id[2:(length(X.id))]))
  #change the observation id
  Xb.id = cumsum(changes.in.X.id)
  #get the different groups of the id
  groups = unique(Xb.id)
  ng = length(groups) #get the number of groups
  T2 = rep(NA, ng) #preallocate memory to the statistic T2
  Z = NULL #preallocate memory for sns

  i = 1 # initialize the group index of the observation id vector
  Yb = Y
  if(!is.null(Yb)){
    Yb = Yb[!is.na(Yb)] # initialize reference sample (remove na values)
  }

  UCL = rep(NA, ng)

  alpha <- chart.par[1]
  nv <- ncol(X)
  switch(chart,
         T2 = {
           if(null.dist == "Chi"){
             ucl <- qchisq(1-alpha,nv) #control limit
           }else if(null.dist=="F"){
             M <- 1
             n <- length(Xb.id) / ng
             if(!is.null(Yb)){#if Yb exists
              m <- nrow(Yb) / ng
              M <- ceiling(m / n)
             }
             ucl <- nv*(M+1)*(n-1)/(M*n-M-nv+1)*qf(1-alpha, nv, M*n-M-nv+1) #control limit
           }
         }
  )

  while (i <= ng) { # repeat until the total groups are analized
    Xb = X[which(Xb.id == groups[i]),] # get the observations to evalute from the positions
    ns = SNS::MNS(X = Xb, Y = Yb, theta = theta, Ftheta = Ftheta, scoring = scoring, alignment = alignment, constant = constant) # calculate the normal score
    Zb = ns$Z
    n = nrow(Xb) #get the number of observation per group

    if (i == 1) { # if is the first batch
      updateSample <- TRUE
      T2[i] = 0 # it does not give any information and is considered the reference sample
    }else{
      mu = apply(Zb, 2, mean) #obtain the mean for each variable in the batch
      # check if the subgroup is in control according to each scheme
      # the reference sample is updated
      updateSample <- FALSE
      switch(chart,
         T2 = {
           T2[i] = n*(mu%*%chol2inv(chol(cor(Z, method = "spearman")))%*%mu) #get the T2 statistic

           if (null.dist == "F"){
             M <- M + 1 #add the subgroup
             ucl <- nv*(M+1)*(n-1)/(M*n-M-nv+1)*qf(1-alpha, nv, M*n-M-nv+1) #control limit
           }
           if (T2[i] <= ucl) updateSample <- TRUE
         }
      )
    }
    UCL[i] = ucl

    if ( (updateSample || is.null(Yb)) && !isFixed){# if the subgroup is in control (updateSample change to TRUE)
      Yb = rbind(Yb, Xb) # add to reference sample the new observations
    }
    Z = rbind(Z, Zb) #update the sns

    i = i + 1 # continue with the next group
  }
  output = list(
    coefficients = list(
      n=n,
      chart = chart
    ),
    X = Xb,
    Z = Z,
    T2 = T2,
    X.id = X.id,
    UCL = UCL
  )

  class(output)="msns" # Class definition

  return(output) # return the sequential normal score
}

#' @import graphics
#' @export
plot.msns <- function(x,...){
  par(mar = c(6,6,4,2))
  T2 = x$T2
  o.id = unique(x$X.id) # original id
  chart = coef(x)$chart
  UCL = x$UCL
  difMaxZ = 0
  difMinZ = 0
  switch(chart,
         T2 = {
           difMaxZ = abs(max(T2) - max(UCL))
           difMinZ = abs(min(T2))
         }
  )

  ymax = max(UCL)
  if(difMaxZ > difMinZ){
    ymax = ymax + difMaxZ
  }else{
    ymax = ymax + difMinZ
  }

  ymin = 0
  if(difMaxZ < difMinZ){
    ymin = ymin - difMaxZ
  }else{
    ymin = ymin - difMinZ
  }

  switch(chart,
         T2 = {
           plot(o.id,T2,type="o",lty=2, lwd=1.5,pch=19,xlab="Batch",ylab=expression(T[SNS]^2),
                ylim=c(ymin, ymax),cex.lab=2.5, cex.axis=1.5, cex=2)
         }
  )

  change = 1
  if(o.id[1] > o.id[length(o.id)]){
    change = -1
  }
  lines(o.id, UCL,lt=4, lwd=3)
}



