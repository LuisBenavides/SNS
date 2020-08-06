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
#' @param omit.id vector. Elements of the vector are the id which are omitted in the analysis.
#' @param auto.omit.alarm logical. Determine if OC signals are added (or not) to reference sample. By default is set to TRUE.
#' @export
#' @examples
#' X = cbind(example91$X1, example91$X2)
#' X.id = example91$X1.id
#' msns = MSNS(X, X.id)
MSNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL, scoring = "Z",
                alignment = "unadjusted", constant = NULL, absolute = FALSE,
                chart="T2", chart.par = c(0.005), null.dist="Chi", isFixed = FALSE,
                omit.id = NULL, auto.omit.alarm = TRUE) {

  if (is.null(theta) != is.null(Ftheta)) { # in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  } else if (nrow(X) != length(X.id)) {
    print("ERROR, observations (X) have different length of the observations id (X.id)")
    return()
  }
  omit.id.found = NULL
  if(!is.null(omit.id)){#check which groups are omitted
    ids = unique(X.id)
    omit.id.found = which(ids %in% omit.id)
    omit.id.missing = omit.id[!(omit.id %in% ids)]
    if(!is.null(omit.id.missing) || length(omit.id.missing) > 0 ){
      cat("WARNING, ids to omit not found:", omit.id.missing, "\n")
    }
    if(auto.omit.alarm){
      auto.omit.alarm = FALSE
      print("WARNING, auto.omit.alarm = FALSE make omit.id = NULL to enable.")
    }
    if(is.null(omit.id.found) || length(omit.id.found) == 0){
      print("WARNING, omitted ids not found, OC signals not added to reference sample (auto.omit.alarm = TRUE).")
    }
  }
  # detect the changes in the observation id vector
  changes.in.X.id = c(1, as.numeric(X.id[1:(length(X.id) - 1)] != X.id[2:(length(X.id))]))
  #change the observation id
  Xb.id = cumsum(changes.in.X.id)
  #get the different groups of the id
  groups = unique(Xb.id)

  ng = length(groups) #get the number of groups
  T2 = rep(NA, ng) #preallocate memory to the statistic T2
  Z = NULL #preallocate memory for sns (IC and OC)
  Zm = NULL #preallocate memory for sns (IC)
  i = 1 # initialize the group index of the observation id vector
  Yb = Y
  if(!is.null(Yb)){
    #Yb = Yb[!is.na(Yb),] # initialize reference sample (remove NA values)

    #get the normal scores of the reference sample
    ns = SNS::MNS(X = Yb, Y = NULL, theta = theta, Ftheta = Ftheta, scoring = scoring, alignment = alignment, constant = constant) # calculate the normal score
    Z = ns$Z
    Zm = Z
  }


  UCL = rep(NA, ng)

  alpha <- chart.par[1]
  nv <- ncol(X)

  switch(chart,
         T2 = {
           if(null.dist == "Chi"){
             ucl <- qchisq(1-alpha,nv) #control limit
           }else if(null.dist == "F"){
             if(FALSE){#check depending on if is known mean or not
               M <- 1
               n <- length(Xb.id) / ng
             }

             m <- 30 # check if it is correct
             # m is the number of observations of the reference sample
             if(!is.null(Yb)){#if there is reference sample
              m <- nrow(Yb)
              if (FALSE){#check depending on if is known mean or not
                m <- mnrow(Yb) / ng
                M <- ceiling(m / n)
              }
             }
             # known mean, unknown covariance matrix
             ucl = ((nv*(m-1))/(m-nv))*qf(1-alpha, nv, m-nv) # control limit

             # Montgomery, normal T2
             # unknown mean, unknown covariance matrix
             #ucl <- nv*(M+1)*(n-1)/(M*n-M-nv+1)*qf(1-alpha, nv, M*n-M-nv+1) #control limit
           }
         }
  )

  while (i <= ng) { # repeat until the total groups are analyzed
    Xb = X[which(Xb.id == groups[i]),] # get the observations to evaluate from the positions
    ns = SNS::MNS(X = Xb, Y = Yb, theta = theta, Ftheta = Ftheta, scoring = scoring, alignment = alignment, constant = constant) # calculate the normal score
    Zb = ns$Z
    n = nrow(Xb) #get the number of observation per group
    if (is.null(Yb)){ # if there is not reference sample
      updateSample <- TRUE

      T2[i] = NA # it does not give any information and is considered the reference sample
    }else{
      #by default there is an update
      updateSample <- TRUE

      corZm = cor(Zm)
      #function to get if the matrix is singular

      isNotSingular <- class(try(chol2inv(chol(corZm)),silent=T))[1] == "matrix"

      if(isNotSingular){#check if the matrix Zm is not is Singular
        mu = apply(Zb, 2, mean) #obtain the mean for each variable in the batch
        m = nrow(Yb) #obtain the number of observations of the reference sample

        # check if the subgroup is in control according to each scheme
        # the reference sample is updated
        updateSample <- FALSE
        switch(chart,
           T2 = {
             T2[i] = n*(mu%*%chol2inv(chol(corZm))%*%mu) #get the T2 statistic

             if (null.dist == "F"){
               #M <- M + 1 #add the subgroup
               ucl = ((nv*(m-1))/(m-nv))*qf(1-alpha, nv, m-nv) #control limit
               #ucl <- nv*(M+1)*(n-1)/(M*n-M-nv+1)*qf(1-alpha, nv, M*n-M-nv+1) #control limit
             }
             if (T2[i] <= ucl){
               updateSample <- TRUE
             }
           }
        )
        UCL[i] = ucl
      }
    }

    if (auto.omit.alarm){
      if ((updateSample || is.null(Yb)) && !isFixed){# if the subgroup is in control (updateSample change to TRUE)
          Yb = rbind(Yb, Xb) # add to reference sample the new observations
          Zm = rbind(Zm, Zb) #update the sns in control
      }
    }else{
      if ((!(i %in% omit.id.found) || is.null(Yb)) && !isFixed){#and if the id of the group is not omitted
        Yb = rbind(Yb, Xb) # add to reference sample the new observations
        Zm = rbind(Zm, Zb) #update the sns in control
      }
    }

    Z = rbind(Z, Zb) #update the sns all in control and out of control

    i = i + 1 # continue with the next group
  }
  output = list(
    coefficients = list(
      n=n,
      chart = chart
    ),
    X = X,
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
  zoom = 0.1
  CEX = 0.75
  switch(chart,
         T2 = {
           difMaxZ = abs(max(T2, na.rm = TRUE) - max(UCL[5:length(UCL)], na.rm = TRUE))
           difMinZ = abs(min(T2, na.rm = TRUE))
         }
  )

  ymax = max(UCL[5:length(UCL)], na.rm = TRUE)
  if(difMaxZ > difMinZ){
    ymax = ymax + difMaxZ
  }else{
    ymax = ymax + difMinZ
  }
  ymax = ymax * (1+zoom)
  ymin = 0
  if(difMaxZ < difMinZ){
    ymin = ymin - difMaxZ
  }else{
    ymin = ymin - difMinZ
  }
  ymin = ymin * (1-zoom)
  switch(chart,
         T2 = {
           plot(o.id,T2,type="o",lty=2, lwd=1.5,pch=19,xlab="Batch",ylab=expression(T[SNS]^2),
                ylim=c(ymin, ymax),cex.lab=CEX*2, cex.axis=CEX*1.25, cex=CEX, col = ifelse(T2 > UCL, "red", "black"))
           lines(o.id, T2, lt=2, lwd=CEX*1.5)
         }
  )

  change = 1
  if(o.id[1] > o.id[length(o.id)]){
    change = -1
  }
  lines(o.id, UCL,lt=4, lwd=CEX*1.5, col = "gray48")
}



