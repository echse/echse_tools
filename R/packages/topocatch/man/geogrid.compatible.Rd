\name{geogrid.compatible}
\alias{geogrid.compatible}
\title{Test two geogrids for spatial compatibility}
\usage{
  geogrid.compatible(grid1, grid2)
}
\arguments{
  \item{grid1}{An object of class geogrid.}

  \item{grid2}{Another object of class geogrid.}
}
\value{
  \code{TRUE} of the two grids are spatially compatible and
  \code{FALSE} otherwise.
}
\description{
  The two geogrids are compared with respect to the grid
  dimensions, the coordinates of the origin, as well as the
  spatial resolution of a grid cell.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

