\name{geogrid}
\alias{geogrid}
\title{Construct a geogrid object from a numeric matrix}
\usage{
  geogrid(m, xllcorner, yllcorner, cellsize, nodata_value)
}
\arguments{
  \item{m}{A numeric matrix.}

  \item{xllcorner}{X-coordinate corresponding to left
  margin of the cell in the lower left corner of the matrix
  \code{m}.}

  \item{yllcorner}{Y-coordinate corresponding to bottom
  margin of the cell in the lower left corner of the matrix
  \code{m}.}

  \item{cellsize}{The size of a quadratic grid cell, i.e.
  cell extent in X and Y direction.}

  \item{nodata_value}{A numeric constant to replace
  \code{NA} values in the matrix \code{m}.}
}
\value{
  An object of class geogrid.
}
\description{
  The functions converts a numeric matrix into a geogrid
  object. The geo-coordinates of the lower left corner and
  the cell size need to be specified. \code{NA} values are
  translated into a special numeric value.
}
\note{
  Infinite values must not appear in the matrix \code{m}
  but \code{NA} values are OK (see argument \code{nodata}).
}
\examples{
g= geogrid(matrix(runif(6*8), ncol=8, nrow=6, byrow=TRUE),
  xllcorner=0, yllcorner=0, cellsize=1, nodata=-9999)
geogrid.plot(g,"Random data")
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

