\name{hour_sunrise}
\alias{hour_sunrise}
\title{Time of sunrise}
\usage{
  hour_sunrise(daynum, latitude)
}
\arguments{
  \item{latitude}{Geographic latitude in (decimal)
  degrees.}

  \item{daynum}{Day number as an integer in range [1,
  366].}
}
\value{
  Time of sunrise after 'darkest' point of the night
  (hours).
}
\description{
  Computes the time of sunrise after the 'darkest' point of
  the night. The function is used in the computation of
  solar radiation from the number of sunshine hours/day.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Bremicker (2000): Das Wasserhaushaltsmodell LARSIM,
  Freiburger Schriften zur Hydrologie, Bd. 11, Eq. 3.23.
  Note that in this source units are wrong (angles are in
  radians rather than degrees).
}

