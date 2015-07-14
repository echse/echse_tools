\name{geogrid.fillGaps}
\alias{geogrid.fillGaps}
\title{Fill gaps in a geogrid}
\usage{
  geogrid.fillGaps(grid, removevalue, maxdist = NULL,
    removeAll = TRUE, silent = TRUE)
}
\arguments{
  \item{grid}{A geogrid object with some raster cells to be
  filled with values from the surrounding.}

  \item{removevalue}{Value indicating the gaps, i.e. the
  value to be substituted in the input grid.}

  \item{maxdist}{The maximum distance allowed in the
  nearest-neighbor search.  If the default of \code{NULL}
  is used, the search distance is set to a value being
  greater than the diagonal of the grid. This guarantees a
  successful search as long as the grid contains some valid
  data.}

  \item{removeAll}{A logical value. If \code{TRUE},
  \code{stop} is called if not all cells with a value of
  \code{removevalue} could be substituted, possibly because
  of a too small value of \code{maxdist}. Defaults to
  \code{TRUE}.}

  \item{silent}{Print status info? (logical)}
}
\value{
  An object of class geogrid.
}
\description{
  The function filles gaps in the matrix of a geogrid by
  the nearest-neighbor method.
}
\examples{
# A sample grid
g= geogrid(m=matrix(c(5,5,5,5,0,5,5,5,5), ncol=3, nrow=3, byrow=TRUE),
  xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
g2= geogrid.fillGaps(g,0,1000,TRUE)
print(g$matrix)
print(g2$matrix)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

