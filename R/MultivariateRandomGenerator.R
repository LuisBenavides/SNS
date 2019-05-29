#' @title Multivariate Random Observations Generetor
#' @description Multivariate Random observations generator selected from several distributions with user defined mean and variance.
#' @param n scalar. Number of observations to be generated.
#' @param nv scalar. Number of variables to be generated.
#' @param dists vector of character string. Distribution of each variable. The length mus be the same as the number of variables. Select from:
#' \itemize{
#'   \item{"Normal": Normal distribution (default).}
#'   \item{"Gamma": Gamma distribution.}
#' }
#' @param mu scalar. Expected value of the desired distribution.
#' @param sigma scalar. Standard deviation of the desired distribution.
#' @param s matrix. Correlation matrix of the variables
#' @param dists.par matrix  For each variable (column), specify
#' \itemize{
#'   \item{\code{par.location}: Location parameter of the desired distribution. Default 0.}
#'   \item{\code{par.scale}: Scale parameter of the desired distribution. Default 1.}
#'   \item{\code{par.shape}: Shape parameter of the desired distribution, Default 1.}
#' }
#' The number of columns must be the same as the number of variables.
#' @param correlation scalar. Corralation between variables.
#' @return A matrix \code{x} with \code{n} observations generated following the selected distribution with its parameters.
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' mgetDist(n=5, nv=2, dists=c("Normal", "Normal"),dists.par= matrix(c(0,1,1,0,1,1), ncol=2))
mgetDist <- function(n, nv, mu = 0, sigma=NULL, correlation=0, s=NULL,  dists = NULL, dists.par = NULL) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package \"MASS\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  mus = rep(mu, nv) # mean

  # correlation matrix
  if (is.null(s)){
    s = diag(nv) #create identity matrix
    pos = which(s != 1) #get position of elements different of 1
    s[pos] = s[pos] + correlation #add the correlation to that positions
  }
  X = MASS::mvrnorm(n, mus, s) # multivariate normal random numbers

  U = pnorm(X) # marginal's cdf

  if (is.null(dists)){
    dists = c("Normal", 0, 1, 1)
    for (j in 1:nv){
      if (j == 1){
        dist.mat = as.matrix(dists)
      }else{
        dist.mat = cbind(dist.mat, dists)
      }
    }
  }else{
      dist.mat = rbind(dists, dists.par)
  }

  X = NULL

  for(j in 1:nv){
    dist = as.character(dist.mat[1,j])
    par.location = as.numeric(dist.mat[2,j])
    par.scale = as.numeric(dist.mat[3,j])
    par.shape = as.numeric(dist.mat[4,j])

    #  inverse CDFs of the distribution
    switch(dist,
       "Normal" = {
         Xv = qnorm(U[,j], mean = par.location, sd = par.scale)
       },
       "Gamma" = {
         Xv = qgamma(U[,j], scale = par.scale, shape=par.shape)
       }
    )
    Xv = as.matrix(Xv, ncol=1)
    if(j == 1){
      X = Xv
    }else{
      X = cbind(X, Xv)
    }
  }
  return(as.matrix(X))
}
