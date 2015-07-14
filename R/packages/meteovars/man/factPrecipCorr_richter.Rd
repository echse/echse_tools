\name{factPrecipCorr_richter}
\alias{factPrecipCorr_richter}
\title{Correction factor for \emph{daily} precipitation after RICHTER}
\usage{
  factPrecipCorr_richter(precip, temper, tempCrit_lower,
    tempCrit_upper, protection_index)
}
\arguments{
  \item{precip}{Vector of precipitation data (mm/day).}

  \item{temper}{Vector of air temperatures corresponding to
  \code{precip} (degree Celsius).}

  \item{tempCrit_lower}{Threshold temperature for mixed
  precip / snow (e.g. -0.5).}

  \item{tempCrit_upper}{Threshold temperature for mixed
  precip / rain (e.g. +0.5).}

  \item{protection_index}{An integer value in range [1,4]
  where 1: open field, 2: slightly protected, 3: moderately
  protected, 4: highly protected. This corresponds to
  horizon angles of 2, 5, 9.5, and 16 degrees,
  respectively.}
}
\value{
  Vector of correction factors to be applied to the
  precipitation data.
}
\description{
  The function computes a factor which can be used to
  correct raw precipitation measurements for the effects of
  wind and evaporation. This approach is applicable to
  \emph{daily} data only. Different correction schemes are
  applied to summer and winter time. The distinction
  between seasons is made based on air temperature.
}
\note{
  The original equation is: \eqn{precipCorr = precip + b *
  precip^e}.  If the correction is expressed as a
  correction factor \eqn{f = precipCorr/precip}, this gives
  \eqn{f = 1 + b * precip^{(e-1)}}.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Richer (1995): Berichte des DWD No. 194, p. 66/67 (paper
  copy).
}

