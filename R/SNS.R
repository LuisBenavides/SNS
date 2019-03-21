#####################################################################################
#
#     SEQUENTIAL NORMAL SCORES
#
#     Authors:  Dr. Victor G. Tercero-Gomez
#               Dr. Luis A. Benavices-Vazquez
#     Affiliation: Tecnologico de Monterrey
#
#     Date: August 23, 2018
#     Version: 1.0
#
#     DESCRIPTION
#
#     Transform a vector X into SNS using initial observations Y if available
#     SNS follow the order of X.
#     Procedure follows Conover et al. (2017)
#
#     Function requires NS, NS1 and NS2 to work.
#
#     X: is a numerical vector.
#     Y: is a numerical vector. Y = NULL if undefined.
#     theta: is a constant
#     Ftheta: is a constant between (0,1)
#
#     COMMENTS
#
#     If ties, average ranks are used.
#     If Y = NULL, normal scores are set relative to X.
#
#     EXAMPLE CONDITIONAL WITH REFERENCE SAMPLE
# 
#     Y = c(10,20,30,40,50,60,70,80,90,100)
#     X = c(30, 35, 45)
#     theta = 40
#     Ftheta = 0.5
#     sample.id = c("a", "b", "c")
#     > SNS(X=X, sample.id=sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#     [1] -0.52440051 -0.31863936  0.08964235
# 
#     EXAMPLE UNCONDITIONAL WITH REFERENCE SAMPLE
#
#     Y = c(10,20,30,40,50,60,70,80,90,100)
#     X = c(30, 35, 45)
#     theta = NULL
#     Ftheta = NULL
#     sample.id = c("a", "b", "c")
#     > SNS(X=X, sample.id=sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#     [1] -0.6045853 -0.3186394  0.0000000
#
#     EXAMPLE CONDITIONAL WITHOUT REFERENCE SAMPLE
#     
#     Y = NULL#c(10,20,30,40,50,60,70,80,90,100)
#     X = c(30, 35, 45)
#     theta = 40
#     Ftheta = 0.5
#     sample.id = c("a", "b", "c")
#     > SNS(X=X, sample.id=sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#     [1] -0.6744898 -0.3186394  0.6744898
#
#     EXAMPLE UNCONDITIONAL WITHOUT REFERENCE SAMPLE
#
#     Y = NULL
#     X = c(30, 35, 45)
#     theta = NULL
#     Ftheta = NULL
#     sample.id = c("a", "b", "c")
#     > SNS(X=X, sample.id=sample.id, Y = Y, theta = theta, Ftheta = Ftheta)
#     [1] 0.0000000 0.6744898 0.9674216
#     
#     References
#     Conover, W. J., Tercero-Gomez, V. G., & Cordero-Franco, A. E. (2017).
#     The sequential normal scores transformation. Sequential Analysis, 36(3), 397-414.


SNS <- function(X, X.id, Y = NULL, theta = NULL, Ftheta = NULL){
  if(is.null(theta) != is.null(Ftheta)){ #in case one is NULL and not the other
    print("ERROR, theta or Ftheta missing")
    return()
  }else if(length(X) != length(X.id)){
    print("ERROR, observations (X) have different length of the observations id (X.id)")
    return()
  }
  #detect the changes in the observation id vector
  changes.in.X.id = c(1,as.numeric(X.id[1:(length(X.id)-1)] != id[2:(length(X.id))]))
  X.id = cumsum(changes.in.X.id) # change the observation id
  #get the different groups of the id
  groups = unique(X.id)
  z = rep(NA,length(X)) #preallocate memory to initialize the SNS
  i = 1    #initialize the group index of the observation id vector
  if(is.null(Y)){ #if there is no reference sample
    Xe = X[which(X.id == 1)] #get the first group
    z[i] = NS(X = Xe, Y = NULL, theta = theta, Ftheta = Ftheta) #calculate the normal score
    Y = Xe #the reference sample is the observations of the first group
    i = i + 1
  }
  if(length(groups) > 1){ #if there are more than one group
    Ye = Y   #initialize the previous information to evaluate
    while(i <= length(groups)){ #repeat until the total groups are analized
      Xe = X[which(sample.id == groups[i])]  #get the observations to evalute from the positions
      z[i] = NS(X = Xe, Y = Ye, theta = theta, Ftheta = Ftheta) #calculate the normal score
      Ye = c(Ye, Xe) #add to reference sample the just evaluated observations
      i = i + 1 #continue with the next group
    }
  }
  return(z) # return the sequential normal score
}


Y = NULL
X = c(30, 35, 45)
theta = NULL
Ftheta = NULL
X.id = c("a", "b")
SNS(X=X, X.id=X.id, Y = Y, theta = theta, Ftheta = Ftheta)