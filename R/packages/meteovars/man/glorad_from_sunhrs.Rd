\name{glorad_from_sunhrs}
\alias{glorad_from_sunhrs}
\title{Estimation of daily-average short-wave radiation from sunshine data}
\usage{
  glorad_from_sunhrs(sunhours, daynum, latitude)
}
\arguments{
  \item{sunhours}{Number of sunshine hours in a day.}

  \item{daynum}{Day number as an integer in range [1,
  366].}

  \item{latitude}{Geographic latitude in (decimal)
  degrees.}
}
\value{
  Radiation in W/m2
}
\description{
  The function estimates daily-average short-wave radiation
  from the number of sunny hours in a day.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Dyck and Peschke (1995): Grundlagen der Hydrologie, Verl.
  f. Bauwesen (see Eq. 30.10 at page 30). The empirical
  coefficients were taken from Bremicker (2000): Das
  Wasserhaushaltsmodell LARSIM, Freiburger Schriften zur
  Hydrologie, Bd. 11 (Eq. 3.21 at page 27).
}

