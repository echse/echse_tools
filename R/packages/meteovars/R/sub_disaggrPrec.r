
#' Disaggregation of sums based on a higher-resolution time series
#'
#' The function takes two vectors as input which usually represent the values
#' of \emph{regular} time series. One of the vectors, \code{L}, is of low resolution
#' and the other one, \code{L} is of higher resolution. Thus, \code{length(L)} is
#' less than \code{length(H)}.
#' This function creates a result vector where the \emph{sum} of the values is taken
#' from \code{L} but the \emph{pattern} in terms of the variability in time is
#' taken from \code{H}.
#'
#' @param L A numeric vector representing the values of the lower-resolution
#'   series whose \emph{sum} is to be retained. The sum is only guaranteed to be
#'   retained exactly if argument \code{setZero} is \code{FALSE}.
#' @param H A numeric vector representing the values of the higher-resolution
#'   series whose \emph{pattern} is to be retained. The lenght of \code{H} must
#'   be a multiple of the length of \code{L}, i.e. \code{length(H) %% length(L)}
#'   must evaluate to zero.
#' @param setZero Logical value to control what happens if, for an element in
#'   \code{L}, the sum of the corresponding elements in \code{H} is zero.
#'   If \code{setZero} is \code{FALSE} (default), the value of the element
#'   in \code{L} is 'equally' distributed over the elements of the result vector.
#'   If \code{setZero} is \code{TRUE}, however, those elements of the result
#'   vector are set to zero. See notes and examples.
#'
#' @return A vector of the same lenght as \code{H} containing the result of the
#'   disaggregation.
#'
#' @note A typical application is the disaggregation of rainfall. Using this
#'   method, it is possible, for example, to transfer the hourly pattern
#'   observed at an automatic rain gage to another observation with daily
#'   resolution made at a nearby rain gage. In that case, \code{setZero=FALSE}
#'   is a reasonable default, because a non-zero observation at the daily rain
#'   gage should be preserved, even if the nearby rain gage with hourly data did
#'   not record any precipitation in any 1-hour interval on that particular day.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' H= c(1,2,3,3,2,1,0,0,0,4,0,4)
#' L= c(2,4,6,8)
#' print(round(disaggSums(L, H, setZero=FALSE),3))
#' print(round(disaggSums(L, H, setZero=TRUE),3))

disaggSums= function(L, H, setZero=FALSE) {
  # check args
  if ((length(H) == 0) || (!is.numeric(H)))
    stop("'H' must be a non-empty numeric vector.")
  if ((length(L) == 0) || (!is.numeric(L)))
    stop("'L' must be a non-empty numeric vector.")
  if ((length(H) %% length(L)) != 0)
    stop("Length of 'H' must be a multiple of length of 'L'.")
  nHperL= length(H)/length(L)
  # process
  # -- append interval index to L, then inflate to the length of H but retain original values
   tmp= data.frame(interval=as.integer(gl(n=length(L), k=nHperL, ordered=TRUE)))
   tmp$v= L[tmp$interval]
   L= tmp
   rm(tmp)
  if (nrow(L) != length(H))
    stop("Failed to inflate 'L'. Result is of bad length.")
  # -- aggregate values in H for the interval levels in L
  H= data.frame(v=H, interval= L$interval)
  agg= tapply(X=H$v, INDEX=H$interval, FUN=sum)
  H$sum= agg[H$interval]
  if (nrow(L) != nrow(H))
    stop("Failed to aggregate 'H'. Result is of bad length.")
  # -- compute weights
  H$weight= H$v / H$sum
  if (setZero) {
    H$weight[!is.finite(H$weight)]= 0
  } else {
    H$weight[!is.finite(H$weight)]= 1/nHperL
  }
  # -- apply weights
  res= L$v * H$weight
  # -- check and return result
  sumL= ifelse(setZero, sum(L$v[which(H$sum != 0)])/nHperL, sum(L$v)/nHperL)
  sumR= sum(res)
  if (sumL != sumR)
    stop(paste("Error in result. Sum over elements in 'L' is ",sumL,
      " but sum over result vector is ",sumR,".",sep=""))
  return(res)
}

