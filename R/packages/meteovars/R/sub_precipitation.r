
################################################################################
#' Precipitation correction factor after SEVRUK
#'
#' The function computes a factor which can be used to correct raw precipitation
#'  measurements for wind effects. Evaporation losses are \emph{not} taken into
#'  accout.
#'
#' @param windsp Vector of wind speed values (m/s).
#' @param temper Vector of air temperatures corresponding to \code{windsp} (degree Celsius).
#' @param sensorHeight_precip Height of rain gage above ground (m).
#' @param sensorHeight_windsp Height of wind sensor above ground (m).
#' @param height_zeroWind Height above ground where wind speed is assumed to
#'   approach zero (m).
#'
#' @return Vector of correction factors to be applied to precipitation data
#'   corresponding to the input vectors \code{windsp} and \code{temper}.
#'
#' @note See the file 'meteo.h' in the ECHSE processes folder for details.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Sevruk, 1989.
#'
#' @export

factPrecipCorr_sevruk= function (windsp, temper, sensorHeight_precip,
  sensorHeight_windsp, height_zeroWind) {
  if (length(windsp) != length(temper)) stop("Length of 'windsp' and 'temper' is different.")
  if (height_zeroWind <= 0.) stop("Argument 'height_zeroWind' must be > 0.")
  if (sensorHeight_windsp <= 0.) stop("Argument 'sensorHeight_windsp' must be > 0.")
  windsp= windsp * log(sensorHeight_precip/height_zeroWind) /
    log(sensorHeight_windsp/height_zeroWind)
  f= rep(NA, length(windsp))
  i= which(temper < -27.) 
  if (length(i) > 0) f[i]= 1. + (0.550 * windsp[i]^1.40)
  i= which((temper >= -27.) & (temper < 8.)) 
  if (length(i) > 0) f[i]= 1. + (0.280 * windsp[i]^1.30)
  i= which((temper >= -8.) & (temper < 0.))
  if (length(i) > 0) f[i]= 1. + (0.150 * windsp[i]^1.18)
  i= which(temper >= 0.)
  if (length(i) > 0) f[i]= 1. + (0.015 * windsp[i]^1.00)
  return( f )
}

################################################################################
#' Correction factor for \emph{daily} precipitation after RICHTER
#'
#' The function computes a factor which can be used to correct raw precipitation
#'  measurements for the effects of wind and evaporation. This approach is
#'  applicable to \emph{daily} data only. Different correction schemes are
#'  applied to summer and winter time. The distinction between seasons is made
#'  based on air temperature.
#'
#' @param precip Vector of precipitation data (mm/day).
#' @param temper Vector of air temperatures corresponding to \code{precip} (degree Celsius).
#' @param tempCrit_lower Threshold temperature for mixed precip / snow (e.g. -0.5).
#' @param tempCrit_upper Threshold temperature for mixed precip / rain (e.g. +0.5).
#' @param protection_index An integer value in range [1,4] where 1: open field,
#'   2: slightly protected, 3: moderately protected, 4: highly protected. This
#'   corresponds to horizon angles of 2, 5, 9.5, and 16 degrees, respectively.
#'
#' @return Vector of correction factors to be applied to the precipitation data.
#'
#' @note The original equation is: \eqn{precipCorr = precip + b * precip^e}.
#'   If the correction is expressed as a correction factor
#'   \eqn{f = precipCorr/precip}, this gives \eqn{f = 1 + b * precip^{(e-1)}}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Richer (1995): Berichte des DWD No. 194, p. 66/67 (paper copy).
#'
#' @export

#
factPrecipCorr_richter= function (precip, temper, tempCrit_lower, tempCrit_upper, protection_index)
{
  if (length(precip) != length(temper)) stop("Length of 'precip' and 'temper' is different.")
  if (tempCrit_lower > tempCrit_upper) stop("'tempCrit_lower' should be <= 'tempCrit_upper'.")
  if (!(protection_index %in% 1:4)) stop("'protection_index' must be one of 1, 2, 3, or 4.")
  # Decide whether to use summer/winter parameters based on temperature
  temp_summer= 7.
  # Parameter set (Richter 1995, p. 67)
  e_rain_summer= 0.38
  e_rain_winter= 0.46
  e_mixed=       0.55
  e_snow=        0.82
  b_rain_summer= c(0.345, 0.310, 0.280, 0.245)
  b_rain_winter= c(0.340, 0.280, 0.240, 0.190)
  b_mixed=       c(0.535, 0.390, 0.305, 0.185)
  b_snow=        c(0.720, 0.510, 0.330, 0.210)
  # Initialize
  f= rep(NA, length(precip))
  # Rain (summer)
  i= which((temper > tempCrit_upper) & (temper > temp_summer))
  if (length(i) > 0) f[i]= 1 + b_rain_summer[protection_index] * precip[i]^(e_rain_summer - 1)
  # Rain (winter)
  i= which((temper > tempCrit_upper) & (temper <= temp_summer))
  if (length(i) > 0) f[i]= 1 + b_rain_winter[protection_index] * precip[i]^(e_rain_winter - 1)
  # Mixed
  i= which((temper <= tempCrit_upper) & (temper > tempCrit_lower))
  if (length(i) > 0) f[i]= 1 + b_mixed[protection_index] * precip[i]^(e_mixed - 1)
  # Snow
  i= which(temper <= tempCrit_lower)
  if (length(i) > 0) f[i]= 1 + b_snow[protection_index] * precip[i]^(e_snow - 1)
  # Correct infinite values due to zero precip
  f[which(f == Inf)]= 1.
  return( f )
}

