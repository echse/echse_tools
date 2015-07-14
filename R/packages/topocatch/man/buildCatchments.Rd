\name{buildCatchments}
\alias{buildCatchments}
\title{Identify sub-basins}
\usage{
  buildCatchments(grid_init, code_conflict, grid_flowdir,
    code_undefined)
}
\arguments{
  \item{grid_init}{A geogrid with integer codes used as
  start points in the process of catchment building. This
  geogrid is typically obtained by converting drainage line
  vectors into a grid using \code{\link{linesToGrid}}.}

  \item{code_conflict}{An integer code marking those cells
  in \code{grid_init} which cannot be unequivocally
  assigned to a catchment. Those cells typically contain
  multiple drainage lines and junctions in particular.}

  \item{grid_flowdir}{A geogrid object holding flow
  direction codes which correspond to the elevation model
  as output by \code{\link{flowdir}}.}

  \item{code_undefined}{An integer code to mark those cells
  in the result grid which cannot be assigned to a
  catchment.}
}
\value{
  An object of class geogrid.
}
\description{
  This function identifies sub-basins based on a grid of
  flow direction codes and an initialization grid. The
  initialization grid is typically obtained by converting
  drainage line vectors into a raster.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

