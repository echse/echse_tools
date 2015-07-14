\name{p_baro}
\alias{p_baro}
\title{Estimation of air pressure based on the barometric formula}
\usage{
  p_baro(z, z0 = 0, p0 = 1013.25, t0 = 15, dtdz = 0.00649)
}
\arguments{
  \item{z}{Elevation of the target station where the
  pressure is to be estimated (m).}

  \item{z0}{Elevation of a reference station where air
  pressure and temperature are known (m).}

  \item{p0}{Air pressure at the reference station (hPa).}

  \item{t0}{Air temperature at the reference station
  (degree C).}

  \item{dtdz}{Temperature lapse rate (K/m).}
}
\value{
  Air pressure at the target location.
}
\description{
  This version of the barometric formula allows for
  estimation of the air pressure at a location using the
  pressure at a reference point and information on the
  temperature lapse rate. A linear temperature gradient is
  assumed to exist in the considered atmospheric layer.
}
\note{
  The default values for \code{z0}, \code{p0}, \code{t0} as
  well as for the temperature lapse rate \code{dtdz} refer
  to a convention known as the international standard
  athmosphere.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  \url{http://de.wikipedia.org/wiki/Barometrische_H\%C3\%B6henformel}
}

