\name{solrad_extra}
\alias{solrad_extra}
\title{Solar radiation at the atmosphere's upper limit}
\usage{
  solrad_extra(daynum, latitude)
}
\arguments{
  \item{daynum}{Day number as an integer in range [1,
  366].}

  \item{latitude}{Geographic latitude in (decimal)
  degrees.}
}
\value{
  Radiation in W/m2
}
\description{
  Computes the solar radiation at the atmosphere's upper
  limit. The function is used in the computation of solar
  radiation from the number of sunshine hours/day. The
  results were checked against values reported in Dyck and
  Peschke (1995): Grundlagen der Hydrologie, Verl. f.
  Bauwesen (see page 31).
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Bremicker (2000): Das Wasserhaushaltsmodell LARSIM,
  Freiburger Schriften zur Hydrologie, Bd. 11. See Eq.
  3.22; Note that in this source, units are wrong (angles
  are in radians rather than in degrees). Also, the result
  unit is Wh/m2 (not W/m2). Therefore, the result must be
  divided by 24 to yield W/m2. The value of the solar
  constant was taken from Wikipedia as 1367 W/m2. This
  agrees with the value of 5022 kJ/m2/h, given in Jokiel,
  C. (1995): Gewaesserguetesimulation natuerlicher
  Fliessgewaesser, ISBN 3883458937, 9783883458939..
}

