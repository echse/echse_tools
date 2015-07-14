
################################################################################
#' Declination of the sun
#' 
#' Computes the sun's declination from the day number. The function is used in
#' the computation of solar radiation from the number of sunshine hours/day.
#'
#' @param daynum Day number as an integer in range [1, 366].
#'
#' @return Declination (radians).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Jokiel, C. (1995): Gewaesserguetesimulation natuerlicher Fliessgewaesser,
#'   ISBN 3883458937, 9783883458939. Also used by Bremicker (2000): Das
#'   Wasserhaushaltsmodell LARSIM, Freiburger Schriften zur Hydrologie, Bd. 11
#'   (see Eq.s 3.22 and 3.23).
#'
#' @export

declination= function(daynum) {
  PI= 3.1415
  return( 23.45/180.*PI * cos(2.*PI*(173-daynum)/365.) )
}

################################################################################
#' Time of sunrise
#' 
#' Computes the time of sunrise after the 'darkest' point of the night. The
#' function is used in the computation of solar radiation from the number of
#' sunshine hours/day.
#'
#' @param latitude Geographic latitude in (decimal) degrees.
#' @inheritParams declination
#'
#' @return Time of sunrise after 'darkest' point of the night (hours).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Bremicker (2000): Das Wasserhaushaltsmodell LARSIM, Freiburger
#'   Schriften zur Hydrologie, Bd. 11, Eq. 3.23. Note that in this source units
#'   are wrong (angles are in radians rather than degrees).
#'
#' @export

hour_sunrise= function (daynum, latitude) {
  PI= 3.1415
  tandec= tan(declination(daynum))
  cosdec= cos(declination(daynum))
  tanlat= tan(latitude*PI/180)
  coslat= cos(latitude*PI/180)
  return( 12/PI*acos(tandec*tanlat+0.0145/(cosdec*coslat)) )
}

################################################################################
#' Number of hours between sunrise and sunset
#' 
#' Computes the max. possible number of sunshine hours in a day
#' (hours between sunrise and sunset). The function is used in the computation
#' of solar radiation from the number of sunshine hours/day.
#'
#' @inheritParams hour_sunrise
#'
#' @return Number of hours between sunrise and sunset (hours)
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Bremicker (2000): Das Wasserhaushaltsmodell LARSIM, Freiburger
#'   Schriften zur Hydrologie, Bd. 11, see above Eq. 3.23.
#'
#' @export

max_sunhours = function (daynum, latitude) {
  24 - 2*hour_sunrise(daynum, latitude)
}

################################################################################
#' Solar radiation at the atmosphere's upper limit
#' 
#' Computes the solar radiation at the atmosphere's upper limit. The function is
#' used in the computation of solar radiation from the number of sunshine hours/day.
#' The results were checked against values reported in Dyck and Peschke (1995):
#' Grundlagen der Hydrologie, Verl. f. Bauwesen (see page 31).
#'
#' @inheritParams hour_sunrise
#'
#' @return Radiation in W/m2
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Bremicker (2000): Das Wasserhaushaltsmodell LARSIM, Freiburger
#'   Schriften zur Hydrologie, Bd. 11. See Eq. 3.22; Note that in this source,
#'   units are wrong (angles are in radians rather than in degrees). Also, the
#'   result unit is Wh/m2 (not W/m2). Therefore, the result must be divided by
#'   24 to yield W/m2. The value of the solar constant was taken from Wikipedia
#'   as 1367 W/m2. This agrees with the value of 5022 kJ/m2/h, given in
#'   Jokiel, C. (1995): Gewaesserguetesimulation natuerlicher Fliessgewaesser,
#'   ISBN 3883458937, 9783883458939..
#'
#' @export

solrad_extra = function (daynum, latitude) {
  PI= 3.1415
  sindec= sin(declination(daynum))
  cosdec= cos(declination(daynum))
  sinlat= sin(latitude*PI/180)
  coslat= cos(latitude*PI/180)
  t1= hour_sunrise(daynum, latitude)
  t2= 24 - t1
  return( 1/24 * 1367. * ((t2-t1) * sindec * sinlat +
    12/PI * cosdec * coslat * (sin(PI*t1/12)-sin(PI*t2/12))) )
}

################################################################################
#' Estimation of daily-average short-wave radiation from sunshine data
#' 
#' The function estimates daily-average short-wave radiation from the number of
#' sunny hours in a day.
#'
#' @param sunhours Number of sunshine hours in a day.
#' @inheritParams hour_sunrise
#'
#' @return Radiation in W/m2
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @references Dyck and Peschke (1995): Grundlagen der Hydrologie,
#'   Verl. f. Bauwesen (see Eq. 30.10 at page 30). The empirical coefficients
#'   were taken from Bremicker (2000): Das Wasserhaushaltsmodell LARSIM,
#'   Freiburger Schriften zur Hydrologie, Bd. 11 (Eq. 3.21 at page 27).
#'
#' @export

glorad_from_sunhrs = function (sunhours, daynum, latitude) {
  a= 0.24
  if ((daynum < 90) || (daynum > (305))) {
    b= 0.5
  } else {
    b= 0.55
  }
  return ( solrad_extra(daynum, latitude) *
    (a + b * sunhours/max_sunhours(daynum, latitude)) )
}

