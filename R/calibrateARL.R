calibrateControlLimit <- function(targetARL=NULL, targetMRL=NULL, n, m, theta=NULL, Ftheta=NULL, dist, mu, sigma, dist.par=c(0, 1, 1), initial.par, replicates=50000, chart, progress=FALSE, parallel=FALSE, maxIter=20){
  #Check for errors
  if(is.null(targetARL) && is.null(targetMRL)){
    print("ERROR: Target ARL or target mRL missing")
    return()
  }else if(!is.null(targetARL) && !is.null(targetMRL)){
    print("ERROR: Two targets defined, delete one")
    return()
  }
  p = 0.1
  if(is.null(targetARL)){
    ARL0 = (targetMRL * 1.5)/ 10
  }else{
    ARL0 = targetARL
  }

  switch(chart,
    Shewhart = {
      name.par = "k"
      index.par = 1
    },
    CUSUM = {
      name.par = "h"
      index.par = 2
    },
    EWMA = {
      name.par = "L"
      index.par = 2
    }
  )

  x = rep(NA, 3)
  y = x

  i = 1
  x[i] = initial.par[index.par]
  while(i < maxIter){
    chart.par[index.par] = x[i]
    result = getARLSNS(n=n, m=m, theta=theta, Ftheta=Ftheta, dist=dist, mu=mu, sigma=sigma, dist.par=dist.par, chart=chart, chart.par=chart.par, print.RL = print.RL, replicates=replicates, progress=progress, parallel=parallel, calibrate=TRUE, arl0=targetARL)
    if(!is.null(targetARL)){
      y[i] = result$ARL
      target = targetARL
      name = "ARL"
    }else{
      y[i] = result$MRL
      target = targetMRL
      name = "MRL"
    }

    if(abs(y[i] - target) <= 0.05*target){
      if(progress) cat("Convergence found with", name.par, "=",x[i],"--",name, "=", y[i], "\n", sep=" ")
      output = list(
        objective.function = y[i],
        par.value = x[i],
        found = TRUE
      )
      return(output)
    }else{

      f1 = 0
      f2 = 0
      if(i > 2){
        f1 = x[i] - target
        f2 = x[i-1] - target
      }

      if(f1*f2 < 0){
        x0 = x[i-1]
        x1 = x[i]
        y0 = y[i-1]
        y1 = y[i]
        m = (y1-y0)/(x1-x0)
        b = y0-m*x0
        x2 = (target-b)/m
        x[i+1] = x2
      }else{
        if(y[i] <= target){
          x[i+1] = x[i] * (1 + p)
        }else{
          x[i+1] = x[i] * (1 - p)
        }
        if(progress) cat("obtained=", y[i], " target=", target," Change h=",x[i]," to h=",x[i+1], "\n",sep="")
      }

    }
    i = i + 1
  }

  posMin = which.min(abs(target-y))
  if(progress) cat("Best",name.par,"found ",x[posMin],"--",name, "=", y[posMin], "\n", sep=" ")

  output = list(
    objective.function = y[posMin],
    par.value = x[posMin],
    found = FALSE
  )
  return(output)

}
