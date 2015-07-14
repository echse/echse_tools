\name{max_sunhours}
\alias{max_sunhours}
\title{Number of hours between sunrise and sunset}
\usage{
  max_sunhours(daynum, latitude)
}
\arguments{
  \item{daynum}{Day number as an integer in range [1,
  366].}

  \item{latitude}{Geographic latitude in (decimal)
  degrees.}
}
\value{
  Number of hours between sunrise and sunset (hours)
}
\description{
  Computes the max. possible number of sunshine hours in a
  day (hours between sunrise and sunset). The function is
  used in the computation of solar radiation from the
  number of sunshine hours/day.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Bremicker (2000): Das Wasserhaushaltsmodell LARSIM,
  Freiburger Schriften zur Hydrologie, Bd. 11, see above
  Eq. 3.23.
}

