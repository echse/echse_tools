\name{factPrecipCorr_sevruk}
\alias{factPrecipCorr_sevruk}
\title{Precipitation correction factor after SEVRUK}
\usage{
  factPrecipCorr_sevruk(windsp, temper,
    sensorHeight_precip, sensorHeight_windsp,
    height_zeroWind)
}
\arguments{
  \item{windsp}{Vector of wind speed values (m/s).}

  \item{temper}{Vector of air temperatures corresponding to
  \code{windsp} (degree Celsius).}

  \item{sensorHeight_precip}{Height of rain gage above
  ground (m).}

  \item{sensorHeight_windsp}{Height of wind sensor above
  ground (m).}

  \item{height_zeroWind}{Height above ground where wind
  speed is assumed to approach zero (m).}
}
\value{
  Vector of correction factors to be applied to
  precipitation data corresponding to the input vectors
  \code{windsp} and \code{temper}.
}
\description{
  The function computes a factor which can be used to
  correct raw precipitation measurements for wind effects.
  Evaporation losses are \emph{not} taken into accout.
}
\note{
  See the file 'meteo.h' in the ECHSE processes folder for
  details.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Sevruk, 1989.
}

