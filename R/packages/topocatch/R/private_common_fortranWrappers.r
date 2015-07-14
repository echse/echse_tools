
################################################################################

#' Filling of sinks in a digital elevation model (DEM)
#'
#' Sinks are filled by a non-iterative approach. To make this work, the data are
#' internally rounded to integers after appropriate scaling. For large elevation
#' models the filling of sinks may consume a lot of time (minutes to hours).
#'
#' @param grid The DEM as an object of class geogrid.
#' @param ndigits Number of relevant digits in the raw DEM (integer); see notes
#' @param silent Print status info? (logical)
#'
#' @return The DEM with sinks filled as a list of type grid.
#'
#' @note Filling of sinks is performed on the rounded elevation data after
#'   multiplication with 10^\code{ndigits}. If the original DEM has a resolution
#'   of decimeters, for example, one would use \code{ndigits} = 1 to preserve
#'   that resolution. If an accurracy of 1 meter is tolerable, one would set
#'   \code{ndigits} = 0 in order to speed up the computation.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # Create a simple DEM with some sinks
#' dem= geogrid(m=matrix(c(
#'     3,3,3,3,3,3,-9999,-9999,
#'     3,1,2,2,3,3,3,-9999,
#'     3,3,3,3,3,3,3,3,
#'     5,5,5,5,5,5,2,5,
#'     5,5,2,3,4,5,2,5,
#'     5,5,5,5,5,5,5,5
#'   ), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' dem2= sinkfill(grid=dem, ndigits=0, silent=FALSE)
#' geogrid.plot(dem, "DEM with sinks")
#' geogrid.plot(dem2, "Sinks filled")

sinkfill= function(grid, ndigits, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=ndigits, len=1, type="integer")
  if (ndigits < 0)
    stop(paste("Argument 'ndigits' must be an integer >= 0.",sep=""))
  checkArg(arg=silent, len=1, type="logical")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("sinkfill",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.double(grid$matrix),
    as.double(grid$nodata_value), as.integer(ndigits), as.integer(silent),
    # Outputs
    out=rep(as.double(0), nc*nr))$out
  # Transform result into grid
  outgrid= grid   # Copy header from input
  outgrid$matrix= matrix(output,ncol=nc,nrow=nr,byrow=FALSE)
  return(outgrid)
}

################################################################################

#' Flow directions on a digital elevation model (DEM)
#'
#' For each cell, the flow direction is encoded as an integer in range 1 to 8.
#' Code 1 means North-West and the code increases in clock-wise direction, thus
#' 2= North, 3= North-East, 4= East and so on. If the drainage direction of a
#' cell is undefined, a code of 0 is assigned to that cell. This happens in the
#' case of (1) sinks, (2) missing elevation data in the surrounding of a cell,
#' and (3) flat regions at the margings of the DEM. Note that flat areas in
#' central part of the DEM are handled (by allowing cells to drain into near
#' drained cells).
#'
#' @param grid The \emph{sink-filled} DEM as an object of class geogrid. It is
#'   important that sinks have been filled before, e.g. using
#'   \code{\link{sinkfill}}.
#' @param silent Print status info? (logical)
#'
#' @return A geogrid of flow direction codes (as integers). See notes.
#'
#' @note Undefined values in the input DEM are passed on to the output. The
#'   corresponding nodata value in the output grid is -9999.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' dem= geogrid(matrix(runif(6*8), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata=-9999)
#' fdir= flowdir(dem)
#' geogrid.plot(fdir,"Flow dir. codes")

flowdir=function(grid, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=silent, len=1, type="logical")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("flowdir",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.double(grid$matrix),
    as.double(grid$nodata_value), as.integer(silent),
    # Outputs
    out=rep(as.integer(0), nc*nr))$out
  # Transform result into matrix
  outgrid= grid                # Copy header from input
  outgrid$nodata_value= -9999  # Assigned in the Fortran routine
  outgrid$matrix= matrix(output,ncol=nc,nrow=nr,byrow=FALSE)
  return(outgrid)
}

################################################################################

#' Flow accumulation based on a grid of flow direction codes
#'
#' The flow accumulation is computed as the number of upstream cells (integer).
#' To convert this into units of an area, one has to multiply by the cell area.
#'
#' @param grid A geogrid object with flow direction codes as output by
#'   \code{\link{flowdir}}.
#' @param silent Print status info? (logical)
#'
#' @return A geogrid of flow accumulation values (as integers). See notes.
#'
#' @note If the flow direction codes were derived from a DEM containing sinks,
#'   the function issues an error messages and calls \code{stop}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # Create a simple DEM
#' dem_raw= geogrid(m=matrix(c(
#'     5,5,5,5,5,5,5,5,
#'     5,4,4,4,4,4,4,5,
#'     5,4,3,3,3,3,4,5,
#'     5,4,3,2,2,3,4,5,
#'     5,4,3,2,2,3,4,5,
#'     5,4,3,1,1,3,4,5
#'   ), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' dem= sinkfill(grid=dem_raw, ndigits=0, silent=FALSE)
#' fdir= flowdir(dem)
#' facc= flowacc(fdir)
#' geogrid.plot(facc,"# upstr. cells")

flowacc=function(grid, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=silent, len=1, type="logical")
  # Check range of flow direction codes
  if (!all(grid$matrix %in% c(0:8,grid$nodata_value)))
    stop("Invalid flow direction codes detected in input grid.")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("flowacc",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.integer(grid$matrix),
    as.integer(grid$nodata_value), as.integer(silent),
    # Outputs
    values=rep(as.integer(0), nc*nr),
    errorlevel= as.integer(0))
  # Evaluate error level
  if (output$errorlevel != 0) {
    if (output$errorlevel == 1) {
      stop("Infinite loop due to sinks in the DEM from which flow direction codes were derived.")
    } else {
      stop("Error in Fortran subroutine. Details not available.")
    }
  }
  # Transform result into matrix
  outgrid= grid                # Copy header from input
  outgrid$nodata_value= -9999  # This is save since values are in range 0...Inf
  outgrid$matrix= matrix(output$values,ncol=nc,nrow=nr,byrow=FALSE)
  return(outgrid)
}

################################################################################

#' Concentration time index derived from DEM
#'
#' The concentration time index (CTI) is a simple DEM-derived measure of the
#' travel time of runoff. For a grid cell of the elevation model, the CTI is
#' defined as \eqn{CTI = \sum_{k=1}^{n} X_k / \sqrt{I_k}} where \eqn{X} is the
#' length of a flow path segment, \eqn{i} is the slope of that segment, and
#' \eqn{n} is the total number of segments in the flow path (i.e. the number of
#' cells through which water must flow until it reaches a stream.
#' The CTI of a cell is small if the flow path to is short (small distance to
#' stream) and the average slope along the flow path is high. 
#' The term \eqn{\sqrt{I_k}} originates from Manning's equation. The density of
#' the stream network is controlled by the critical source area (CSA).
#'
#' @param grid_dem A geogrid object representing the \emph{sink-filled} digital
#'   elevation model as output by \code{\link{sinkfill}}.
#' @param grid_flowdir A geogrid object holding flow direction codes which
#'   correspond to the elevation model as output by \code{\link{flowdir}}.
#' @param grid_flowacc A geogrid object with flow accumulation values
#'   corresponding to the elevation model as output by \code{\link{flowacc}}.
#' @param crit_source_area A numeric value defining the critical source area
#'   (CSA). Streams are assumed to form in cells with a catchment whose areal
#'   extent is >= the CSA value. Should be several times larger than the extent
#'   of a single grid cell. 
#' @param dz_min A small numeric value > 0 representing a minimum elevation
#'   difference between adjacent cells of the elevation model. It prevents the
#'   slope and thus the CTI to become infinite in flat areas. Defaults to 0.1.
#' @param silent Print status info? (logical)
#'
#' @return A geogrid of CTI values.
#'
#' @note Undefined CTI values (due to undefined cells in the DEM or the grid of
#'   flow direction) are set to -9999.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # Create a simple DEM
#' dem_raw= geogrid(m=matrix(c(
#'     5,5,5,5,5,5,5,5,
#'     5,4,4,4,4,4,4,5,
#'     5,4,3,3,3,3,4,5,
#'     5,4,3,2,2,3,4,5,
#'     5,4,3,2,2,3,4,5,
#'     5,4,3,1,1,3,4,5
#'   ), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' dem= sinkfill(grid=dem_raw, ndigits=0, silent=FALSE)
#' fdir= flowdir(grid=dem)
#' facc= flowacc(grid=fdir)
#' cti= concTimeIndex(dem, fdir, facc, crit_source_area=2, dz_min=0.1)
#' geogrid.plot(cti, "CT index")

concTimeIndex= function(grid_dem, grid_flowdir, grid_flowacc, crit_source_area,
  dz_min=0.1, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid_flowdir))
    stop("Argument 'grid_flowdir' is not a valid geogrid.")
  if (!is.geogrid(grid_dem))
    stop("Argument 'grid_dem' is not a valid geogrid.")
  if (!is.geogrid(grid_flowacc))
    stop("Argument 'grid_flowacc' is not a valid geogrid.")
  checkArg(arg=crit_source_area, len=1, type="numeric")
  checkArg(arg=dz_min, len=1, type="numeric")
  checkArg(arg=silent, len=1, type="logical")
  # Check grid headers for compatibility  
  if ((!geogrid.compatible(grid_dem, grid_flowdir)) ||
      (!geogrid.compatible(grid_dem, grid_flowacc)))
    stop("Input grids not spatially compatible.")
  # Get dimensions
  nc= grid_dem$ncols
  nr= grid_dem$nrows
  # Call Fortran
  output= .Fortran("conctimeindex",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr),
    as.integer(grid_flowdir$matrix), as.integer(grid_flowdir$nodata_value),
    as.double(grid_dem$matrix), as.double(grid_dem$nodata_value),
    as.integer(grid_flowacc$matrix),
    as.double(grid_dem$cellsize), as.double(crit_source_area), as.double(dz_min),
    as.integer(silent),
    # Outputs
    matrix_cti=rep(as.double(0.), nc*nr), errorlevel=as.integer(0))
  # Evaluate error level
  if (output$errorlevel != 0) {
    if (output$errorlevel == 1) {
      stop("Value of 'crit_source_area' too small. Must be >= cell extent.")
    } else if (output$errorlevel == 2) {
      stop("Value of 'dz_min' must be positive.")
    } else if (output$errorlevel == 3) {
      stop("Missing elevation data for cell(s) with a valid flow direction code.")
    } else {
      stop("Error in Fortran subroutine. Details not available.")
    }
  } else {
    # Transform result into matrix
    outgrid= grid_dem            # Copy header from input
    outgrid$nodata_value= -9999. # Assigned in the Fortran routine
    outgrid$matrix= matrix(output$matrix_cti,ncol=nc,nrow=nr,byrow=FALSE)
    if (!is.geogrid(outgrid)) stop("Result is not a valid grid.")
    return(outgrid)
  }
}

################################################################################

#' Compute flow path vectors
#'
#' Flow path vectors are identified from a grid of the flow accumulation which
#' is derived from the flow direction grid.
#'
#' @param grid_flowdir A geogrid object holding flow direction codes computed
#'   from an elevation model using \code{\link{flowdir}}.
#' @param grid_flowacc A geogrid object with flow accumulation values
#'   corresponding to the same elevation model as output by \code{\link{flowacc}}.
#' @param crit_source_area A numeric value defining the critical source area
#'   (CSA). Streams are assumed to form in cells with a catchment whose areal
#'   extent is >= the CSA value. Should be several times larger than the extent
#'   of a single grid cell. 
#' @param x_inBasin A x-coordinate which is known to be located inside the
#'   river basin of interest. If \code{NULL}, the center point of the input
#'   grid is used.
#' @param y_inBasin A y-coordinate which is known to be located inside the
#'   river basin of interest. If \code{NULL}, the center point of the input
#'   grid is used.
#' @param silent Print status info? (logical)
#'
#' @return A list with two components 'shpTable' and 'id_outlet' The first
#'   component is a data frame with three columns 'id', 'x', and 'y'. The second
#'   component is an integer value, representing the id of the flow path at the
#'   system's downstream end.
#'
#' @note Flow path are returned for a single river basin only, which is defined
#'   by \code{x_inBasin} and \code{y_inBasin}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # A sample DEM
#' dem_raw= geogrid(m=matrix(c(
#'     5,5,5,5,5,5,5,5,
#'     5,4,4,4,4,4,4,5,
#'     5,4,3,3,3,3,4,5,
#'     5,4,3,2,2,3,4,5,
#'     5,4,3,2,1,3,4,5,
#'     5,4,3,2,1,3,4,5
#'   ), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' dem= sinkfill(grid=dem_raw, ndigits=0, silent=FALSE)
#' fdir= flowdir(grid=dem)
#' facc= flowacc(grid=fdir)
#' streams= flowPaths(grid_flowdir=fdir, grid_flowacc=facc, crit_source_area=2,
#'   x_inBasin=4.5, y_inBasin=1.5)
#' streams.plot= function() {
#'   for (id in unique(streams$shpTable$id)) {
#'     i= which(streams$shpTable$id == id)
#'     lines(x=streams$shpTable$x[i], y=streams$shpTable$y[i],
#'       lty=1+2*(id == streams$id_outlet))
#'   }
#' }
#' geogrid.plot(dem, "DEM + paths", fun=streams.plot)

flowPaths= function(grid_flowdir, grid_flowacc, crit_source_area,
  x_inBasin=NULL, y_inBasin=NULL, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid_flowdir))
    stop("Argument 'grid_flowdir' is not a valid geogrid.")
  if (!is.geogrid(grid_flowacc))
    stop("Argument 'grid_flowacc' is not a valid geogrid.")
  checkArg(arg=crit_source_area, len=1, type="numeric")
  if (!is.null(x_inBasin))
    checkArg(arg=x_inBasin, len=1, type="numeric")
  if (!is.null(y_inBasin))
    checkArg(arg=y_inBasin, len=1, type="numeric")
  checkArg(arg=silent, len=1, type="logical")
  # Check compatibility of grid headers
  if (!geogrid.compatible(grid_flowdir, grid_flowacc))
    stop("Input grids not spatially compatible.")
  # Get dimensions
  nc= grid_flowdir$ncols
  nr= grid_flowdir$nrows
  # Set in-basin location if not given
  if (is.null(x_inBasin)) x_inBasin= grid_flowdir$xllcorner + grid_flowdir$cellsize * grid_flowdir$ncols/2
  if (is.null(y_inBasin)) y_inBasin= grid_flowdir$yllcorner + grid_flowdir$cellsize * grid_flowdir$nrows/2
  # Call Fortran
  output= .Fortran("flowpaths",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.integer(grid_flowdir$matrix), as.integer(grid_flowacc$matrix),
    as.double(grid_flowdir$xllcorner), as.double(grid_flowdir$yllcorner), as.double(grid_flowdir$cellsize),
    as.double(crit_source_area), as.double(x_inBasin), as.double(y_inBasin),
    as.integer(silent),
    # Outputs
    ncoords=as.integer(0),
    path_id= rep(as.integer(0), nc*nr*2),
    path_x= rep(as.double(0.), nc*nr*2),
    path_y= rep(as.double(0.), nc*nr*2),
    id_outlet=as.integer(0),
    errorlevel=as.integer(0))
  # Evaluate error level
  if (output$errorlevel != 0) {
    if (output$errorlevel == 1) {
      stop("Value of 'crit_source_area' too small. Must be >= cell extent.")
    } else if (output$errorlevel == 2) {
      stop("Bad in-basin coordinates. Location is outside the grid's extent.")
    } else if (output$errorlevel == 3) {
      stop("Bad in-basin coordinates. Location has non-usable flow direction code.")
    } else {
      stop("Error in Fortran subroutine. Details not available.")
    }
  }
  # Collect results in a single table of the actual length
  if (output$ncoords < 2)
    stop("Result table contains less than 2 records.")
  tab= data.frame(
    id= output$path_id[1:output$ncoords],
    x= output$path_x[1:output$ncoords],
    y= output$path_y[1:output$ncoords]
  )
  # Sort by id
  tab= tab[sort.list(tab$id),]
  # Check computed outlet id
  if (!(output$id_outlet %in% tab$id))
    stop(paste("Computed ID of system outlet (",output$id_outlet,") is invalid.",sep=""))
  # Return list
  return(list(shpTable= tab, id_outlet= output$id_outlet))
}

################################################################################

#' Vector-to-raster conversion for lines
#'
#' Converts a set of lines into a raster.
#'
#' @param ncol Number of columns in the output raster
#' @param nrow Number of rows in the output raster
#' @param xllcorner X-coordinate of the output raster's Western margin
#' @param yllcorner Y-coordinate of the output raster's Southern margin
#' @param cellsize The size of a quadratic grid cell, i.e. cell extent in X and
#'   Y direction.
#' @param code_undefined An integer code to be used as a background value for
#'   the output raster. This code must not be present in the 'id' column of the
#'   \code{lineTable} argument. It is used as the 'nodata' value of the result
#'   grid.
#' @param code_conflict Integer code assigned to those cells of the output
#'   raster touched by multiple lines. This code must not be present in the 'id'
#'   column of the \code{lineTable} argument.
#' @param nbuffer An integer value to control the thickness of the rasterized
#'   lines. The default value of 0 produces lines being just 1 cell wide. Values
#'   of \code{nbuffer}=1, 2, 3 result in line withs of 3, 5, 7 cells,
#'   respectively.
#' @param lineTable A data frame holding the lines' IDs and coordinates. There
#'   must be three columns named 'id', 'x', and 'y'. The 'id' column must be of
#'   type integer while the values in 'x' and 'y' are expected to be numeric.
#' @param silent Print status info? (logical)
#'
#' @return An object of class geogrid.
#'
#' @note There is no guarantee that all lines being present in \code{lineTable}
#'   can be transferred into the output raster. This is especially so if a
#'   line's length is in the order of the raster's cell size or smaller. In
#'   those cases, the function reports the ID(s) of the problematic line(s) and
#'   calls \code{stop}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' arcs= data.frame(
#'   id= c(1,1,1,2,2,2,3,3,3,3,3,3),
#'   x=  c(1,2,3,5,4,3,3,2,2,3,4,5),
#'   y=  c(5,4,3,5,4,3,3,2,1,1,1,1)
#' )
#' raster= linesToGrid(ncol=5, nrow=5, xllcorner=0.5, yllcorner=0.5, cellsize=1,
#'   lineTable=arcs)
#' lines.plot= function() {
#'   for (id in unique(arcs$id)) {
#'     i= which(arcs$id == id)
#'     lines(x=arcs$x[i], y=arcs$y[i])
#'   }
#' }
#' geogrid.plot(raster,"Gridded lines",fun=lines.plot)

linesToGrid= function(ncol, nrow, xllcorner, yllcorner, cellsize, lineTable,
  code_undefined=-9999, code_conflict=-1, nbuffer=0, silent=TRUE) {
  # Check matrix constructor data
  checkArg(arg=ncol, len=1, type="integer")
  checkArg(arg=nrow, len=1, type="integer")
  checkArg(arg=xllcorner, len=1, type="numeric")
  checkArg(arg=yllcorner, len=1, type="numeric")
  checkArg(arg=cellsize, len=1, type="numeric")
  checkArg(arg=code_undefined, len=1, type="numeric")
  checkArg(arg=code_conflict, len=1, type="numeric")
  checkArg(arg=nbuffer, len=1, type="integer")
  checkArg(arg=lineTable, len=NULL, type="data.frame")
  if (!nrow(lineTable))
    stop("Table in argument 'lineTable' must not have zero rows.")
  cols= c("x","y","id")
  if (!all(cols %in% names(lineTable)))
    stop(paste("Argument 'lineTable' must be a data frame with columns '",paste(cols,collapse="', '"),"'.",sep=""))
  if (!all(is.numeric(unlist(lineTable[,]))))
    stop("Argument 'lineTable' must be a data frame with only numeric data.")
  if (any(lineTable$id == code_undefined))
    stop("Value of argument 'code_undefined' is already used as a line ID.")
  if (any(lineTable$id == code_conflict))
    stop("Value of argument 'code_conflict' is already used as a line ID.")
  checkArg(arg=silent, len=1, type="logical")
  # Call Fortran
  output= .Fortran("lines2raster",
    # Inputs
    as.integer(ncol), as.integer(nrow),
    as.double(xllcorner), as.double(yllcorner), as.double(cellsize),
    as.integer(code_undefined), as.integer(code_conflict), as.integer(nbuffer),
    as.integer(silent),
    as.integer(nrow(lineTable)), as.double(lineTable$x), as.double(lineTable$y),
    as.integer(lineTable$id),
    # Outputs
    out=rep(as.integer(0), ncol*nrow))$out
  # Stop if not all lines do appear in the grid (usually because
  # they are too short and the grid resolution is too coarse)
  lineIDs= unique(lineTable$id)
  disappeared= which(!(lineIDs %in% unique(output)))
  if (length(disappeared) > 0) {
    stop(paste(length(disappeared)," of ",
      length(lineIDs)," lines do not appear in the output grid. These lines",
      " may be too short for the grid's resolution. The affected line IDs follow: '",
      paste(lineIDs[disappeared],collapse="', '"),"'.",sep=""))
  }
  # Transform result into matrix
  outgrid= geogrid(m=matrix(output,ncol=ncol,nrow=nrow,byrow=FALSE),
    xllcorner=xllcorner, yllcorner=yllcorner, cellsize=cellsize,
    nodata_value=code_undefined)
  if (!is.geogrid(outgrid)) stop("Result is not a valid geogrid.")
  return(outgrid)
}

################################################################################

#' Identify sub-basins
#'
#' This function identifies sub-basins based on a grid of flow direction codes
#' and an initialization grid. The initialization grid is typically obtained by
#' converting drainage line vectors into a raster.
#'
#' @param grid_init A geogrid with integer codes used as start points in the
#'   process of catchment building. This geogrid is typically obtained by
#'   converting drainage line vectors into a grid using
#'   \code{\link{linesToGrid}}.
#' @param code_conflict An integer code marking those cells in \code{grid_init}
#'   which cannot be unequivocally assigned to a catchment. Those cells 
#'   typically contain multiple drainage lines and junctions in particular.
#' @param grid_flowdir A geogrid object holding flow direction codes which
#'   correspond to the elevation model as output by \code{\link{flowdir}}.
#' @param code_undefined An integer code to mark those cells in the result grid
#'   which cannot be assigned to a catchment.
#'
#' @return An object of class geogrid.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch

buildCatchments= function(grid_init, code_conflict, grid_flowdir, code_undefined) {
  # Check args
  if (!is.geogrid(grid_init)) stop("Argument 'grid_init' is not a valid geogrid.")
  checkArg(arg=code_conflict, len=1, type="numeric")
  if (!is.geogrid(grid_flowdir)) stop("Argument 'grid_flowdir' is not a valid geogrid.")
  checkArg(arg=code_undefined, len=1, type="numeric")
  # Check input grids for compatibility
  if (!geogrid.compatible(grid_init, grid_flowdir))
    stop("Input grids are not spatially compatible.")
  # Get dimensions
  nc= grid_init$ncols
  nr= grid_init$nrows
  # Call Fortran
  output= .Fortran("buildcatchments",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr),
    as.integer(grid_init$nodata_value), as.integer(code_conflict), as.integer(grid_init$matrix),
    as.integer(grid_flowdir$nodata_value), as.integer(grid_flowdir$matrix),
    as.integer(code_undefined),
    # Outputs
    out=rep(as.integer(0), nc*nr))$out
  # Transform result into matrix
  outgrid= grid_init           # Copy header from input
  outgrid$matrix= matrix(output,ncol=nc,nrow=nr,byrow=FALSE)
  if (!is.geogrid(outgrid)) stop("Result is not a valid geogrid.")
  return(outgrid)
}

################################################################################

#' Fill gaps in a geogrid
#'
#' The function filles gaps in the matrix of a geogrid by the nearest-neighbor
#' method.
#'
#' @param grid A geogrid object with some raster cells to be filled with values
#'   from the surrounding.
#' @param removevalue Value indicating the gaps, i.e. the value to be
#'   substituted in the input grid.
#' @param maxdist The maximum distance allowed in the nearest-neighbor search.
#'   If the default of \code{NULL} is used, the search distance is set to a
#'   value being greater than the diagonal of the grid. This guarantees a
#'   successful search as long as the grid contains some valid data. 
#' @param removeAll A logical value. If \code{TRUE}, \code{stop} is called if
#'   not all cells with a value of \code{removevalue} could be substituted,
#'   possibly because of a too small value of \code{maxdist}. Defaults to
#'   \code{TRUE}.
#' @param silent Print status info? (logical)
#'
#' @return An object of class geogrid.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # A sample grid
#' g= geogrid(m=matrix(c(5,5,5,5,0,5,5,5,5), ncol=3, nrow=3, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' g2= geogrid.fillGaps(g,0,1000,TRUE)
#' print(g$matrix)
#' print(g2$matrix)

geogrid.fillGaps= function(grid, removevalue, maxdist=NULL, removeAll=TRUE, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid.")
  checkArg(arg=removevalue, len=1, type="numeric")
  if (is.null(maxdist)) {
    maxdist= grid$ncols*grid$cellsize + grid$nrows*grid$cellsize + grid$cellsize
  } else {
    checkArg(arg=maxdist, len=1, type="numeric")
  }
  checkArg(arg=removeAll, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("nnfill",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.double(grid$cellsize),
    as.integer(grid$nodata_value), as.integer(grid$matrix),
    as.integer(removevalue), as.double(maxdist), as.integer(silent),
    # Outputs
    out=rep(as.integer(0), nc*nr))$out
  # Check for non-removed values
  nleft= sum(output == removevalue)
  if ((nleft > 0) && removeAll) stop(paste(nleft," cells with value ",
    removevalue," remained unchanged.",sep=""))
  # Transform result into matrix
  outgrid= grid           # Copy header from input
  outgrid$matrix= matrix(output,ncol=nc,nrow=nr,byrow=FALSE)
  if (!is.geogrid(outgrid)) stop("Result is not a valid geogrid.")
  return(outgrid)
}

################################################################################

#' Extract values from a geogrid
#'
#' The function returns the values of a geogrid's matrix corresponding to a set
#' of locations.
#'
#' @param grid A geogrid object.
#' @param xy_table A data frame defining the sampling points. There must be at
#'   least two columns named 'x' and 'y'.
#' @param silent Print status info? (logical)
#'
#' @return A numeric vector whose length is the number rows in \code{xy_table}.
#'   For sampling points being located outside the grid's extension, the grid's
#'   \emph{nodata} value (i.e. \code{grid$nodata_value}) is returned.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # A sample grid
#' g= geogrid(m=matrix(1:9, ncol=3, nrow=3, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' p= data.frame(x= c(0.5,2.5), y= c(0.5,2.5))
#' print(geogrid.valuesAtPoints(g,p))

geogrid.valuesAtPoints= function(grid, xy_table, silent=TRUE) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid.")
  checkArg(arg=xy_table, len=NULL, type="data.frame")
  if (!nrow(xy_table)) stop("Table in argument 'xy_table' must not have zero rows.")
  cols= c("x","y")
  if (!all(cols %in% names(xy_table)))
    stop(paste("Argument 'xy_table' must be a data frame with columns '",paste(cols,collapse="', '"),"'.",sep=""))
  if (!all(is.numeric(unlist(xy_table[,]))))
    stop("Argument 'xy_table' must be a data frame with only numeric data.")
  checkArg(arg=silent, len=1, type="logical")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("gridvaluesatpoints",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(nc), as.integer(nr), as.double(grid$xllcorner), as.double(grid$yllcorner),
    as.double(grid$cellsize), as.double(grid$nodata_value), as.double(grid$matrix),
    as.integer(nrow(xy_table)), as.double(xy_table$x), as.double(xy_table$y),
    as.integer(silent),
    # Outputs
    out=as.double(rep(0, nrow(xy_table))))$out
  return(output)
}

################################################################################

#' Create geogrid holding the index of the nearest points in a set of points
#'
#' The function returns a geogrid where the value of each cell is the index
#' of the nearest point in a list of points.
#'
#' @param ncol Number of columns in the output raster
#' @param nrow Number of rows in the output raster
#' @param xllcorner X-coordinate of the output raster's Western margin
#' @param yllcorner Y-coordinate of the output raster's Southern margin
#' @param cellsize The size of a quadratic grid cell, i.e. cell extent in X and
#'   Y direction.
#' @param xy_table A data frame defining the sampling points. There must be at
#'   least two columns named 'x' and 'y'.
#' @param silent Print status info? (logical)
#'
#' @return A geogrid of the specified dimensions filled with integer cell values
#'   between 1 and \code{nrow(xy_table)}. The nodata value is set to zero.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' points= data.frame(id=1:4, x=c(2,5,2,5), y=c(3,3,8,8))
#' g= geogrid.indexOfNearestPoints(ncol=6, nrow=10, xllcorner=0, yllcorner=0,
#'   cellsize=1, xy_table=points)
#' print(g$matrix)

geogrid.indexOfNearestPoints= function(ncol, nrow, xllcorner, yllcorner,
  cellsize, xy_table, silent=TRUE) {
  # Check args
  checkArg(arg=ncol, len=1, type="integer")
  checkArg(arg=nrow, len=1, type="integer")
  checkArg(arg=xllcorner, len=1, type="numeric")
  checkArg(arg=yllcorner, len=1, type="numeric")
  checkArg(arg=cellsize, len=1, type="numeric")
  checkArg(arg=xy_table, len=NULL, type="data.frame")
  if (!nrow(xy_table)) stop("Table in argument 'xy_table' must not have zero rows.")
  cols= c("x","y")
  if (!all(cols %in% names(xy_table)))
    stop(paste("Argument 'xy_table' must be a data frame with columns '",paste(cols,collapse="', '"),"'.",sep=""))
  if (!all(is.numeric(unlist(xy_table[,]))))
    stop("Argument 'xy_table' must be a data frame with only numeric data.")
  checkArg(arg=silent, len=1, type="logical")
  # Call Fortran
  output= .Fortran("gridnearestpoints",
    # Inputs (matrix is passed as vector of concatenated COLUMNS (not rows))
    as.integer(ncol), as.integer(nrow), as.double(xllcorner), as.double(yllcorner),
    as.double(cellsize),
    as.integer(nrow(xy_table)), as.double(xy_table$x), as.double(xy_table$y),
    as.integer(silent),
    # Outputs
    out=rep(as.integer(0), ncol*nrow))$out
  # Transform result into matrix
  outgrid= geogrid(m=matrix(output,ncol=ncol,nrow=nrow,byrow=FALSE),
    xllcorner=xllcorner, yllcorner=yllcorner, cellsize=cellsize, nodata_value=0)
  if (!is.geogrid(outgrid)) stop("Result is not a valid geogrid.")
  return(outgrid)
}

################################################################################

#' Modify/reclassify values of a geogrid object
#'
#' This function allows for re-classification of the values in a geogrid's
#' matrix.
#'
#' @param grid A geogrid object.
#' @param lower A numeric vector defining the lower limits of the classes in the
#'   result grid.
#' @param upper A numeric vector defining the upper limits of the classes in the
#'   result grid. Must be of the same length as \code{lower}.
#' @param new A numeric vector holding the values to be assigned to the classes.
#'   Must be of the same length as \code{lower} and \code{upper}.
#' @param includeLower Should values being equal those in \code{lower} become
#'   part of the classes? (logical)
#' @param includeUpper Should values being equal those in \code{upper} become
#'   part of the classes? (logical)
#'
#' @return An object of class geogrid.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib topocatch
#'
#' @examples
#' # A sample grid
#' g= geogrid(m=matrix(1:9, ncol=3, nrow=3, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata_value=-9999)
#' print(geogrid.reclass(g, lower=c(1,5), upper=c(5,9), new=c(0,1)))

geogrid.reclass= function(grid, lower, upper, new=0.5*(lower+upper), includeLower=TRUE, includeUpper=TRUE) {
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=lower, len=NULL, type="numeric")
  checkArg(arg=upper, len=NULL, type="numeric")
  checkArg(arg=new, len=NULL, type="numeric")
  if ((length(lower) != length(upper)) || (length(lower) != length(new)))
    stop("Input vectors 'lower', 'upper', and 'new' must be of the same length.")
  checkArg(arg=includeLower, len=1, type="logical")
  checkArg(arg=includeUpper, len=1, type="logical")
  # Get dimensions
  nc= grid$ncols
  nr= grid$nrows
  # Call Fortran
  output= .Fortran("vectreclass",
    # Inputs
    as.integer(nc*nr), as.double(grid$matrix), as.integer(length(new)),
    as.double(lower), as.double(upper), as.double(new),
    as.integer(includeLower), as.integer(includeUpper),
    # Output
    out=rep(as.double(0), nc*nr))$out
  # Transform result into matrix
  outgrid= grid           # Copy header from input
  outgrid$matrix= matrix(output,ncol=nc,nrow=nr,byrow=FALSE)
  if (!is.geogrid(outgrid)) stop("Result is not a valid geogrid.")
  return(outgrid)
}


