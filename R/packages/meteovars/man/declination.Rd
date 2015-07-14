\name{declination}
\alias{declination}
\title{Declination of the sun}
\usage{
  declination(daynum)
}
\arguments{
  \item{daynum}{Day number as an integer in range [1,
  366].}
}
\value{
  Declination (radians).
}
\description{
  Computes the sun's declination from the day number. The
  function is used in the computation of solar radiation
  from the number of sunshine hours/day.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Jokiel, C. (1995): Gewaesserguetesimulation natuerlicher
  Fliessgewaesser, ISBN 3883458937, 9783883458939. Also
  used by Bremicker (2000): Das Wasserhaushaltsmodell
  LARSIM, Freiburger Schriften zur Hydrologie, Bd. 11 (see
  Eq.s 3.22 and 3.23).
}

