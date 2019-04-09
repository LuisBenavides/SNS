#' Alignment of the data
#' @description Align the monitoring sample \code{X} and the reference sample \code{Y}.
#' @param X vector. Monitoring sample.
#' @param Y vector. Reference sample.
#' @param alignment character string. Aligment of the data \code{X} and \code{Y}. Select from
#' \itemize{
#'   \item "unadjusted": nothing is sustracte from \code{X} and \code{Y} (default).
#'   \item "overallmean": overall mean is sustracted from \code{X} and \code{Y}.
#'   \item "overallmedian": overall median is sustracted from \code{X} and \code{Y}.
#'   \item "samplemean": mean from corresponding group (\code{X} and \code{Y}) is sustracted from its corresponing vector.
#'   \item "samplemedian": median from corresponding group (\code{X} and \code{Y}) is sustracted from its corresponing vector.
#'   \item "referencemean": mean from \code{Y} is sustracted from \code{X} and \code{Y}.
#'   \item "referencemedian": median from \code{Y} is sustracted from \code{X} and \code{Y}.
#'   \item "constantvalue": a constant value is sustracted from \code{X} and \code{Y}.
#' }
#' @param constant scalar. Only used when the \code{alignment} is selected "constantvalue". Default \code{NULL}.
#' @param absolute logical. If \code{TRUE}, the absolute aligned values are obtained. (Default \code{FALSE})
#' @export
#' @examples
#' X = c(30, 45, 50)
#' Y = c(20, 22, 25, 30, 70)
#' dataAlignment(X,Y)
#'
dataAlignment <- function(X, Y,
                          alignment="unadjusted", constant=NULL, absolute=FALSE){
  #Alignment
  switch (alignment,
    unadjusted = {
      x.adjusted = X
      y.adjusted = Y
    },
    overallmean = {
      omean = mean(c(Y,X))
      x.adjusted = X - omean
      y.adjusted = Y - omean
    },
    overallmedian = {
      omedian = median(c(Y, X))
      x.adjusted = X - omedian
      y.adjusted = Y - omedian
    },
    samplemean = {
      x.adjusted = X - mean(X)
      y.adjusted = Y - mean(Y)
    },
    samplemedian = {
      x.adjusted = X - median(X)
      y.adjusted = Y - median(Y)
    },
    referencemean = {
      rmean = mean(Y)
      x.adjusted = X - rmean
      y.adjusted = Y - rmean
    },
    referencemedian = {
      rmedian = median(Y)
      x.adjusted = X - rmedian
      y.adjusted = Y - rmedian
    },
    constantvalue={
      if(is.null(constant)){
        print("ERROR: Must specify constant value.")
        return()
      }
      x.adjusted = X - constant
      y.adjusted = Y - constant
    })

  if(length(y.adjusted) == 0){
    y.adjusted = NULL
  }

  if(absolute){
    x.adjusted = abs(x.adjusted)
    y.adjusted = abs(y.adjusted)
  }

  output = list(
    X = x.adjusted,
    Y = y.adjusted
  )
  return(output)
}
