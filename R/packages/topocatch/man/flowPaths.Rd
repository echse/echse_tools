\name{flowPaths}
\alias{flowPaths}
\title{Compute flow path vectors}
\usage{
  flowPaths(grid_flowdir, grid_flowacc, crit_source_area,
    x_inBasin = NULL, y_inBasin = NULL, silent = TRUE)
}
\arguments{
  \item{grid_flowdir}{A geogrid object holding flow
  direction codes computed from an elevation model using
  \code{\link{flowdir}}.}

  \item{grid_flowacc}{A geogrid object with flow
  accumulation values corresponding to the same elevation
  model as output by \code{\link{flowacc}}.}

  \item{crit_source_area}{A numeric value defining the
  critical source area (CSA). Streams are assumed to form
  in cells with a catchment whose areal extent is >= the
  CSA value. Should be several times larger than the extent
  of a single grid cell.}

  \item{x_inBasin}{A x-coordinate which is known to be
  located inside the river basin of interest. If
  \code{NULL}, the center point of the input grid is used.}

  \item{y_inBasin}{A y-coordinate which is known to be
  located inside the river basin of interest. If
  \code{NULL}, the center point of the input grid is used.}

  \item{silent}{Print status info? (logical)}
}
\value{
  A list with two components 'shpTable' and 'id_outlet' The
  first component is a data frame with three columns 'id',
  'x', and 'y'. The second component is an integer value,
  representing the id of the flow path at the system's
  downstream end.
}
\description{
  Flow path vectors are identified from a grid of the flow
  accumulation which is derived from the flow direction
  grid.
}
\note{
  Flow path are returned for a single river basin only,
  which is defined by \code{x_inBasin} and
  \code{y_inBasin}.
}
\examples{
# A sample DEM
dem_raw= geogrid(m=matrix(c(
    5,5,5,5,5,5,5,5,
    5,4,4,4,4,4,4,5,
    5,4,3,3,3,3,4,5,
    5,4,3,2,2,3,4,5,
    5,4,3,2,1,3,4,5,
    5,4,3,2,1,3,4,5
  ), ncol=8, nrow=6, byrow=TRUE),
  xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
dem= sinkfill(grid=dem_raw, ndigits=0, silent=FALSE)
fdir= flowdir(grid=dem)
facc= flowacc(grid=fdir)
streams= flowPaths(grid_flowdir=fdir, grid_flowacc=facc, crit_source_area=2,
  x_inBasin=4.5, y_inBasin=1.5)
streams.plot= function() {
  for (id in unique(streams$shpTable$id)) {
    i= which(streams$shpTable$id == id)
    lines(x=streams$shpTable$x[i], y=streams$shpTable$y[i],
      lty=1+2*(id == streams$id_outlet))
  }
}
geogrid.plot(dem, "DEM + paths", fun=streams.plot)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

