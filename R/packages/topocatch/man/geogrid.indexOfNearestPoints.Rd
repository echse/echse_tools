\name{geogrid.indexOfNearestPoints}
\alias{geogrid.indexOfNearestPoints}
\title{Create geogrid holding the index of the nearest points in a set of points}
\usage{
  geogrid.indexOfNearestPoints(ncol, nrow, xllcorner,
    yllcorner, cellsize, xy_table, silent = TRUE)
}
\arguments{
  \item{ncol}{Number of columns in the output raster}

  \item{nrow}{Number of rows in the output raster}

  \item{xllcorner}{X-coordinate of the output raster's
  Western margin}

  \item{yllcorner}{Y-coordinate of the output raster's
  Southern margin}

  \item{cellsize}{The size of a quadratic grid cell, i.e.
  cell extent in X and Y direction.}

  \item{xy_table}{A data frame defining the sampling
  points. There must be at least two columns named 'x' and
  'y'.}

  \item{silent}{Print status info? (logical)}
}
\value{
  A geogrid of the specified dimensions filled with integer
  cell values between 1 and \code{nrow(xy_table)}. The
  nodata value is set to zero.
}
\description{
  The function returns a geogrid where the value of each
  cell is the index of the nearest point in a list of
  points.
}
\examples{
points= data.frame(id=1:4, x=c(2,5,2,5), y=c(3,3,8,8))
g= geogrid.indexOfNearestPoints(ncol=6, nrow=10, xllcorner=0, yllcorner=0,
  cellsize=1, xy_table=points)
print(g$matrix)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

