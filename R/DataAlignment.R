#' Alignment of the data
#'
#' Align the monitoring sample X and the reference sample Y
#'
#' @param X is a numerical vector.
#' @param Y is a numerical vector.
#' @param alignment is the aligment of the data
#' unadjusted:
#' overallmean:
#' overallmedian:
#' samplemean:
#' samplemedian:
#' referencemean:
#' referencemedian:
#' constantvalue:
#' @param constant is a numeric value
#' @export
#' @examples
#' X = c(30, 45, 50)
#' Y = c(20, 22, 25, 30, 70)
#' dataAlignment(X,Y)
#'
dataAlignment <- function(X, Y, alignment="unadjusted", constant=NULL){
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
  output = list(
    X = x.adjusted,
    Y = y.adjusted
  )
  return(output)
}
