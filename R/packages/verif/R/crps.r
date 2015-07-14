
#' CRPS of a forecast (ensemble) of a continuous variable
#'
#' The function computes the continous ranked probability score (CRPS) for an
#' ensemble (or deterministic) forecast of a continous variable.
#'
#' @param obs Vector of observations of length m (m being the number of cases).
#' @param fct Matrix holding predicted values corresponding to \code{obs}. There
#'   must be as many rows as there are elements in \code{obs}. For a deterministic
#'   forecast, the matrix has a single column only. In the case of an ensemble
#'   forecast, each column represents a member.
#'
#' @return For compatibility with the 'crps' function from the verification
#'   package, a list with two components is returned:
#'   \item{crps}{The computed crps for each single case (vector).}
#'   \item{CRPS}{The crps averaged over all cases (scalar).}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Hersbach, H. (2000): Decomposition of the Continuous Ranked Probability
#'   Score for Ensemble Prediction Systems; Weather and forecasting, Vol. 15, pages
#'   559-570 (used Eqn.s 5, 24, 25, 26, 27, and info from text at page 563)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' Compare with the function 'crps' from package 'verification'. This function
#' doesn't use the values of the individual members but only mean and standard
#' deviation of the forecasts.
#'
#' Generate observations
#' ncases= 1000
#' obs= runif(ncases)
#' # Generate mean and sdevs of enemble forecasts
#' fct_means= runif(ncases)
#' fct_sdevs= runif(ncases)
#' # Generate ensemble forecasts with choosen mean and sdev
#' nmembs= 20
#' fct= matrix(NA, nrow=ncases, ncol=nmembs)
#' for (i in 1:ncases) {
#'   # Draw a perfect (non-random) sample from the normal distribution
#'   fct[i,]= qnorm(p=seq(from=0.01, to=0.99, length.out=nmembs),
#'     mean=fct_means[i], sd=fct_sdevs[i])
#' }
#' # Compute crps using the discrete member values
#' ver1= crps_hersbach(obs, fct)
#' # Compute crps using only mean and sdev
#' require(verification)
#' ver2= crps(obs, cbind(fct_means, fct_sdevs))
#' # Compare
#' plot(ver1$crps, ver2$crps, col="grey")
#' points(ver1$CRPS, ver2$CRPS, pch=20, col="black")
#' abline(0,1)
#'}

crps_hersbach = function (obs, fct) {

  if (!is.vector(obs)) stop("'obs' must be a vector.")
  ncases= length(obs)
  # If fct is a vector (deterministic fcst) make it a single-column matrix
  if (!(is.vector(fct) || is.matrix(fct))) stop("'fct' must be a matrix or vector.")
  if ((!is.matrix(fct)) && (length(fct) == ncases)) {
    fct= as.matrix(fct, nrow= ncases, byrow=T)
  }
  if (nrow(fct) != length(obs)) stop("Number of rows in 'fct' <> number of elements in 'obs'.")
  nmembs= ncol(fct)

  # Check input
  if (any(is.na(obs))) stop("NA values in 'obs' detected.")
  if (any(is.na(fct))) stop("NA values in 'fct' detected.")

  # Sort members (smallest member = 1st column, largest member = last column)
  for (icase in 1:ncases) {
    fct[icase,]= sort(fct[icase,])
  }

  # Compute the crps for each case
  p= c(0, (1:nmembs)/nmembs)            # probabilities for the bins (Eqn. 23)
  crps_case= rep(0., ncases)            # initialization
  for (icase in 1:ncases) {
    for (i in 1:(nmembs+1)) {
      # Set limits of the bins (based on member values extended to -Inf/Inf)
      # See text below Eqn. 23
      if (i == 1) {                 # limits of the bin below smallest member value
        xmin= -Inf
        xmax= fct[icase,1]
      } else if (i == (nmembs+1)) { # limits of the bin above largest member value
        xmin= fct[icase,nmembs]
        xmax= Inf
      } else {                      # limits of intermediate bins
        xmin= fct[icase,i-1]
        xmax= fct[icase,i]
      }
      # Variables a and b for non-outliers (Eqn. 26)
      if (obs[icase] > xmax) {
        if (i == 1) {               # first bin: only outliers would contribute to crps!
          a= 0                      # See text above Eqn. 27
          b= 0
        } else {
          a= xmax - xmin
          b= 0
        }
      } else if (obs[icase] < xmin) {
        if (i == (nmembs+1)) {      # last bin: only outliers would contribute to crps!
          a= 0                      # See text above Eqn. 27
          b= 0
        } else {
          a= 0
          b= xmax - xmin
        }
      } else {
        a= obs[icase] - xmin
        b= xmax - obs[icase]
      }
      # Variables a and b for outliers (Eqn. 27)
      if ((i == 1) && (obs[icase] <= xmax)) {
        a= 0
        b= xmax - obs[icase]
      }
      if ((i == (nmembs+1)) && (obs[icase] >= xmin)) {
        a= obs[icase] - xmin
        b= 0
      }
      # crps for the current case (Eqn. 24 & 25)
      crps_case[icase]= crps_case[icase] + a * p[i]^2 + b * (1-p[i])^2
    }
  }
  # Return crps averaged over the cases
  return(list(
    crps= crps_case,
    CRPS= sum(crps_case) / ncases,   # Eqn. 5
    cases= ncases,
    members= nmembs
  ))
}


