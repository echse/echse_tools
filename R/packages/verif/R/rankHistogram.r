
#' Rank histogram for an ensemble forecast
#'
#' The rank histrogram is a convenient measure of verification applicable
#' to ensemble forcasts. It is a measure of how well the ensemble spread
#' represents the true variability (uncertainty) of the observations.
#'
#' @param obs A numeric vector of observations.
#' @param fct A matrix representing the ensemble forecasts to be verified. There
#'   must be as many rows in \code{fct} as there are elements in vector \code{obs}.
#'   Each column of the matrix represents an ensemble member.
#' @param plot If \code{TRUE}, the rank histogram is plotted.
#'
#' @return A data frame with the columns 'bin' and 'count'. The number of bins
#'   equals the number of ensemble members + 1. The 'zero bin' counts the number
#'   of observed values being less than the smallest ensemble member.
#'
#' @note The rank histrogram may be interpreted as follows
#'   (\url{http://www.cawcr.gov.au/projects/verification}):
#'
#'  Flat - ensemble spread about right to represent forecast uncertainty
#'
#'  U-shaped - ensemble spread too small, many observations falling outside the
#'  extremes of the ensemble
#'
#'  Dome-shaped - ensemble spread too large, most observations falling near the
#'  center of the ensemble
#'
#'  Asymmetric - ensemble contains bias
#'
#'  Attention: A flat rank histogram does not necessarily indicate a good forecast, it only
#'  measures whether the observed probability distribution is well represented
#'  by the ensemble.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references \url{http://www.cawcr.gov.au/projects/verification};
#'   Hamill, T.M., 2001: Interpretation of rank histograms for verifying ensemble forecasts. Mon. Wea. Rev., 129, 550-560.
#'
#' @export
#'
#' @examples
#' # Forecast and observations draw from same distribution (the perfect case)
#' nevents=5000
#' nmembers=20
#' fct= matrix(runif(nevents*nmembers), ncol=nmembers)
#' obs= runif(nevents)
#' rankHistogram(obs=obs,fct=fct,plot=TRUE)

rankHistogram= function(
  obs,
  fct,
  plot=FALSE
) {
  # Check args
  tryCatch({
    stopifnot(is.numeric(obs),length(obs)>0,!any(is.na(obs)))
    stopifnot(is.matrix(fct),is.numeric(fct),nrow(fct)==length(obs),!any(is.na(fct)))
  }, error= function(e) {
    stop(paste("Error in arguments. Details:",e))
  })
  # Algorithm as described at http://www.cawcr.gov.au/projects/verification
  # 1. At every observation (or analysis) point rank the N ensemble members from
  #    lowest to highest. This represents N+1 possible bins that the observation
  #    could fit into, including the two extremes
  # 2. Identify which bin the observation falls into at each point
  # 3. Tally over many observations to create a histogram of rank.
  rhist=data.frame(bin=0:ncol(fct), count=0)
  for (i in 1:length(obs)) {
    f= sort(fct[i,])
    bin= approx(x=f, y=1:length(f), xout=obs[i], method="constant", f=0,
      ties="ordered", yleft=0, yright=ncol(fct))$y
    rhist$count[bin+1]= rhist$count[bin+1] + 1
  }
  if (plot) {
    barplot(height=rhist$count/sum(rhist$count), names.arg=rhist$bin,
      xlab="Bin", ylab="Relative frequency")
  }
  return(rhist)
}


