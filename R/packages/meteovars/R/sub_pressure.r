################################################################################

#' Estimation of air pressure based on the barometric formula
#'
#' This version of the barometric formula allows for estimation of the air
#' pressure at a location using the pressure at a reference point and
#' information on the temperature lapse rate. A linear temperature gradient
#' is assumed to exist in the considered atmospheric layer.
#'
#' @param z Elevation of the target station where the pressure is to be estimated (m).
#' @param z0 Elevation of a reference station where air pressure and temperature
#'   are known (m).
#' @param p0 Air pressure at the reference station (hPa). 
#' @param t0 Air temperature at the reference station (degree C).
#' @param dtdz Temperature lapse rate (K/m).
#'
#' @return Air pressure at the target location.
#'
#' @note The default values for \code{z0}, \code{p0}, \code{t0} as well as for
#'   the temperature lapse rate \code{dtdz} refer to a convention known as the
#'   international standard athmosphere.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references \url{http://de.wikipedia.org/wiki/Barometrische_H\%C3\%B6henformel}
#'
#' @export

p_baro= function(z, z0=0., p0=1013.25, t0=15., dtdz=0.00649)
{
  return( p0 * (1 - dtdz*(z-z0)/(t0+273.16))^(0.02896*9.807/8.314/dtdz) )
}

