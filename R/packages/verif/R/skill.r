
#' Binary forecast verification measures applied to continous variables
#'
#' The function computes the probability of detection (POD), the false alarm
#' rate (FAR), and the equitable threat score (ETS) for a series of observations
#' and one or more predictions.
#'
#' @param x A data frame holding observations and the corresponding
#'   predictions.
#' @param colObs Name of the column in \code{x} holding the observations.
#' @param colsSim A vector of names identifying the columns in \code{x} that
#'   contain the predictions.
#' @param thresholds A numeric vector defining the thresholds for which scores
#'   are computed. The default uses selected quantiles of the observations.
#' @param nMin The minimum number of observed (POD) or predicted (FAR) events to
#'   accept the results. If, for a particular threshold value, the number of
#'   events is less than \code{nMin}, the result for this threshold is set to
#'   \code{NA}. ETS is set to \code{NA} if either of the number of observed and
#'   predicted events is less than \code{nMin}.
#'
#' @return A list with the components 'POD', 'FAR', and 'ETS'. Each of theses is
#'   a data frame with a column 'threshold' followed by as many columns as
#'   present in \code{colsSim}. Column names are copied from that argument.
#'
#' @note For POD and FAR, the range is 0 (worst) to 1 (perfect). The range of
#'    ETS is -1/3 (worst) to 1 (perfect). An ETS of 0 indicates no skill.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#'
#' @references \url{http://www.cawcr.gov.au/projects/verification}
#'
#' @export
#'
#' @examples
#' d= data.frame(obs=runif(1000), sim1=runif(1000), sim2=runif(1000))
#' skillAsBinary(x=d, colsSim=c("sim1","sim2"),  thresholds=c(0.5,0.75,0.9))

skillAsBinary= function(
  x,
  colObs="obs",
  colsSim=c("sim"),
  thresholds=quantile(x[,colObs], probs=c(0.9, 0.95, 0.99)),
  nMin=5
) {
  # Check args
  stopifnot(is.data.frame(x))
  stopifnot(colObs %in% names(x))
  stopifnot(all(colsSim %in% names(x)))
  if (colObs %in% colsSim) stop("Bad selection of column name(s).")
  # Initialize result tables
  pod= data.frame(threshold= thresholds, matrix(NA, ncol=length(colsSim),
    nrow=length(thresholds)))
  names(pod)[2:ncol(pod)]= colsSim
  far= data.frame(threshold= thresholds, matrix(NA, ncol=length(colsSim),
    nrow=length(thresholds)))
  names(far)= names(pod)
  ets= data.frame(threshold= thresholds, matrix(NA, ncol=length(colsSim),
    nrow=length(thresholds)))
  names(ets)= names(pod)
  # Compute
  for (sim in colsSim) {
    for (i in 1:length(thresholds)) {
      hits=    sum((x[,colObs] >= thresholds[i]) & (x[,sim] >= thresholds[i]))
      misses=  sum((x[,colObs] >= thresholds[i]) & (x[,sim] < thresholds[i]))
      falarms= sum((x[,colObs] < thresholds[i])  & (x[,sim] >= thresholds[i]))
      # POD and FAR
      pod[i,sim]= ifelse((hits + misses) >= nMin, hits / (hits + misses), NA)
      far[i,sim]= ifelse((hits + falarms) >= nMin, falarms / (hits + falarms), NA)
      # ETS
      hits_random= (hits+misses)*(hits+falarms)/nrow(x)
      ets[i,sim]= ifelse((((hits + misses) >= nMin) && ((hits + falarms) >= nMin)),
        (hits - hits_random) / (hits - hits_random + misses + falarms), NA)
    }
  }
  return(list(POD=pod, FAR=far, ETS=ets))
}



