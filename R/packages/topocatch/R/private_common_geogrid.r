
################################################################################

#' Test whether an object is of class geogrid
#'
#' A valid geogrid object is a of class geogrid and consists of a list with
#' seven named components 'ncols', 'nrows', 'xllcorner', 'yllcorner',
#' 'cellsize', 'nodata_value', and 'matrix'.
#'
#' @param x The object to be tested.
#'
#' @return A logical value.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

is.geogrid= function(x) {
  if (is.null(class(x))) return(FALSE)
  if (class(x)[1] != "geogrid") return(FALSE)
  if (!is.list(x)) return(FALSE)
  if (!identical(names(x), c("ncols","nrows","xllcorner","yllcorner","cellsize",
    "nodata_value","matrix"))) {
    warning("Missing or unexpected list components detected.")
    return(FALSE)
  }
  if (!all(is.numeric(c(x$ncols,x$nrows,x$xllcorner,x$yllcorner,
    x$cellsize,x$nodata_value)))) {
    warning("Non-numeric data in header info.")
    return(FALSE)
  }
  if (!all(sapply(list(x$ncols,x$nrows,x$xllcorner,x$yllcorner,
    x$cellsize,x$nodata_value),length) == 1)) {
    warning("Non-scalar data in header info.")
    return(FALSE)
  }
  if (!is.matrix(x$matrix)) {
    warning("Matrix component is not a matrix.")
    return(FALSE)
  }
  if (ncol(x$matrix) != x$ncol) {
    warning("Number of columns in matrix doesn't match with header info.")
    return(FALSE)
  }
  if (nrow(x$matrix) != x$nrow) {
    warning("Number of rows in matrix doesn't match with header info.")
    return(FALSE)
  }
  if ((x$nrow * x$ncol) == 0) {
    warning("Matrix is empty.")
    return(FALSE)
  }
  if (!is.numeric(x$matrix)) {
    warning("Matrix must consist of numeric values.")
    return(FALSE)
  }
  if (!all(is.finite(x$matrix))) {
    warning("Matrix must not contain values like NA, Inf, ect.")
    return(FALSE)
  }
  return(TRUE)
}

################################################################################

#' Test two geogrids for spatial compatibility
#'
#' The two geogrids are compared with respect to the grid dimensions, the
#' coordinates of the origin, as well as the spatial resolution of a grid cell.
#'
#' @param grid1 An object of class geogrid.
#' @param grid2 Another object of class geogrid.
#'
#' @return \code{TRUE} of the two grids are spatially compatible and
#'   \code{FALSE} otherwise.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

geogrid.compatible= function(grid1, grid2) {
  if (!is.geogrid(grid1)) stop("Input 'grid1' is not a valid geogrid.")
  if (!is.geogrid(grid2)) stop("Input 'grid2' is not a valid geogrid.")
  if (grid1$ncols != grid2$ncols) return(FALSE)
  if (grid1$nrows != grid2$nrows) return(FALSE)
  if (grid1$cellsize != grid2$cellsize) return(FALSE)
  if (abs(grid1$xllcorner - grid2$xllcorner) > grid1$cellsize/1.e06) return(FALSE)
  if (abs(grid1$yllcorner - grid2$yllcorner) > grid1$cellsize/1.e06) return(FALSE)
  return(TRUE)
}

################################################################################

#' Read a geogrid object from a file
#'
#' The expected file format is ESRI's ASCII grid format, i.e. there must be
#' six header lines with keywords 'ncols', 'nrows', 'xllcorner', 'yllcorner',
#' 'cellsize', and 'nodata_value'. The matrix of values starts at line seven.
#'
#' @param file Name/path of the ASCII grid file.
#'
#' @return An object of class geogrid.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Create a sample file
#' g= geogrid(matrix(runif(6*8), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata=-9999)
#' f= tempfile()
#' g= geogrid.writeAscii(grid=g, file=f)
#' # Read and display
#' g= geogrid.readAscii(f)
#' geogrid.plot(g,"Random data")
#' file.remove(f)

geogrid.readAscii= function(file) {
  # Check args
  checkFileIn(file)
  # Read header
  h= read.table(file=file, header=FALSE, sep="", nrows=6, colClasses=c("character","numeric"))
  if (ncol(h) != 2) stop(paste("Cannot read header of grid file '",file,"'.",sep=""))
  names(h)= c("item","value")
  h$item= tolower(h$item)
  items= c("ncols","nrows","xllcorner","yllcorner","cellsize","nodata_value")
  if (!all(items %in% h$item)) stop(paste("Missing header key(s) in grid file '",file,"'.",sep=""))
  header= h$value[match(items,h$item)]
  names(header)= items
  # Read matrix  
  m= scan(file=file, what=double(0), nmax=header["ncols"]*header["nrows"], sep="",
    skip=6, nlines=header["nrows"], na.strings=NA, quiet=T)
  # Return
  result=vector("list", 7)
  result[1:6]= header
  names(result)[1:6]= items
  result[[7]]= matrix(m, ncol=header["ncols"], nrow=header["nrows"], byrow=T)
  names(result)[7]= "matrix"
  class(result)= "geogrid"
  return(result)
}

################################################################################

#' Construct a geogrid object from a numeric matrix
#'
#' The functions converts a numeric matrix into a geogrid object. The
#' geo-coordinates of the lower left corner and the cell size need to be
#' specified. \code{NA} values are translated into a special numeric value.
#'
#' @param m A numeric matrix.
#' @param xllcorner X-coordinate corresponding to left margin of the cell in the
#'   lower left corner of the matrix \code{m}.
#' @param yllcorner Y-coordinate corresponding to bottom margin of the cell in
#'   the lower left corner of the matrix \code{m}.
#' @param cellsize The size of a quadratic grid cell, i.e. cell extent in X and
#'   Y direction.
#' @param nodata_value A numeric constant to replace \code{NA} values in the
#'   matrix \code{m}.
#'
#' @return An object of class geogrid.
#'
#' @note Infinite values must not appear in the matrix \code{m} but \code{NA}
#'   values are OK (see argument \code{nodata}).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' g= geogrid(matrix(runif(6*8), ncol=8, nrow=6, byrow=TRUE),
#'   xllcorner=0, yllcorner=0, cellsize=1, nodata=-9999)
#' geogrid.plot(g,"Random data")

geogrid= function(m, xllcorner, yllcorner, cellsize, nodata_value) {
  # Check args
  if ((!is.matrix(m)) || (!is.numeric(m)))
    stop("Argument 'm' must be a numeric matrix.")
  if (any(is.infinite(m)))
    stop("Argument 'm' must not contain infinite values.")
  checkArg(arg=xllcorner, len=1, type="numeric")
  checkArg(arg=yllcorner, len=1, type="numeric")
  checkArg(arg=cellsize, len=1, type="numeric")
  checkArg(arg=nodata_value, len=1, type="numeric")
  # Set nodata value
  m[is.na(m)]= nodata_value
  # Construct and return geogrid
  g= list(ncols= ncol(m), nrows= nrow(m), xllcorner= xllcorner,
    yllcorner= yllcorner, cellsize= cellsize, nodata_value= nodata_value,
    matrix=m)
  class(g)= "geogrid"
  if (!is.geogrid(g))
    stop("Failed to construct geogrid object from input.")
  return(g)
}


################################################################################

#' Export a geogrid object to a file
#'
#' The grid is exported to a file in ESRI's ASCII grid format, i.e. there are
#' six header lines with keywords 'ncols', 'nrows', 'xllcorner', 'yllcorner',
#' 'cellsize', and 'nodata_value'. The matrix of values starts at line seven.
#'
#' @param grid A geogrid object.
#' @param file Name/path of the result file.
#' @param replace Is is OK to replace an existing output file? Defaults to
#'   \code{FALSE}.
#'
#' @return \code{NULL}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

geogrid.writeAscii= function(grid, file, replace=FALSE) {
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=replace, len=1, type="logical")
  checkFileOut(file, replace)
  h= data.frame(
    names= c("ncols","nrows","xllcorner","yllcorner","cellsize","nodata_value"),
    values= c(grid$ncols,grid$nrows,grid$xllcorner,grid$yllcorner,grid$cellsize,grid$nodata_value))
  write.table(x=h, file=file, sep="\t", col.names=F, row.names=F, quote=F)
  write.table(x=grid$matrix, file=file, sep=" ", col.names=FALSE, row.names=FALSE, quote=FALSE, append=T)
  return(invisible(NULL))
}


################################################################################
# Internal functions to determine the limits and center of a grid cell

cell_xmin= function(col,xll,cellsize) xll + (col-1)*cellsize
cell_xmax= function(col,xll,cellsize) xll + (col)*cellsize
cell_ymin= function(row,nrows,yll,cellsize) yll + (nrows-row)*cellsize
cell_ymax= function(row,nrows,yll,cellsize) yll + (nrows-(row-1))*cellsize
cell_xavg= function(col,xll,cellsize) {xll + (col-1) * cellsize + 0.5 * cellsize}
cell_yavg= function(row,nrows,yll,cellsize) {yll + (nrows-row) * cellsize + 0.5 * cellsize}

################################################################################

#' Plot a geogrid object
#'
#' The geogrid's matrix is converted into a raster image that can be either
#' send to a file or the screen.
#'
#' @param grid A geogrid object.
#' @param title A title string to appear above the plot's legend.
#' @param randomColors If \code{TRUE}, a random color is assigned to each unique
#'   value. Otherwise, a continuous color scale is used.
#' @param len The number of entries in the legend of continuous plots. The value is
#'   only used if \code{randomColors} is \code{TRUE}.
#' @param fun A function containing secondary plot commands (like
#'   \code{points} or \code{lines} or \code{text}) to draw additional
#'   objects or annotations on the grid's surface. Defaults to a function
#'   with an empty body.
#' @param file A file name for the output. If \code{NULL} (default), the output
#'   appears on the screen only.
#' @param type A file extension string.  Currently 'png' or 'jpg'.  If
#'   \code{file} is \code{NULL}, the value of this argument is ignored.
#' @param replace Is is OK to replace an existing output file? If \code{file} is
#'   \code{NULL}, the value of this argument is ignored.
#' @param ... Optional arguments to be passed to \code{fun}.
#'
#' @return \code{NULL}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

geogrid.plot= function(grid, title="", randomColors=FALSE, len=if(randomColors) 0 else 5,
  fun=function(){}, file=NULL, type="png", replace=FALSE, ...
) {
  # Check args
  if (!is.geogrid(grid)) stop("Argument 'grid' is not a valid geogrid object.")
  checkArg(arg=title, len=1, type="character")
  checkArg(arg=randomColors, len=1, type="logical")
  checkArg(arg=len, len=1, type="integer")
  if ((!randomColors) && (len < 2))
    len=2
  checkArg(arg=fun, len=1, type="function")
  types= c("png", "jpg")
  if ((length(type)!=1) || (!(type %in% types)))
    stop(paste("Argument 'type' must be one of '",paste(types,collapse="', '"),"'.",sep=""))
  checkArg(arg=fun, len=1, type="function")
  # Start
  nc= grid$ncols
  nr= grid$nrows
  cz= grid$cellsize
  xll= grid$xllcorner
  yll= grid$yllcorner
  legend_frac= 0.33
  # Basic plot
  op= par(no.readonly=T)
  par(mar=c(3,3,2,0.5))
  plot(x=c(xll, xll+nc*cz + nc*cz*legend_frac), y=c(yll, yll+nr*cz), type="n",
    bty="n", xaxt="n", yaxt="n", xlab="x", ylab="y", asp=1)
  # x-axes
  t= axTicks(side=1)
  t= t[(t >= xll) & (t <= xll+nc*cz)]
  labs= format(t,scientific=FALSE)
  axis(side=1, pos=yll, at=t, labels=labs, lwd=0, lwd.ticks=1)
  axis(side=3, pos=yll+nr*cz, at=t, labels=FALSE, lwd=0, lwd.ticks=1)
  # y-axes
  t= axTicks(side=2)
  t= t[(t >= yll) & (t <= yll+nr*cz)]
  labs= format(t,scientific=FALSE)
  axis(side=2, pos=xll, at=t, labels=labs, lwd=0, lwd.ticks=1)
  axis(side=4, pos=xll+nc*cz, at=t, labels=FALSE, lwd=0, lwd.ticks=1)
  # Map
  if (!randomColors) {
    clr_nodata= "grey"
    # Actual limits
    zmin= min(grid$matrix[grid$matrix != grid$nodata_value])
    zmax= max(grid$matrix[grid$matrix != grid$nodata_value])
    leg= pretty(c(zmin,zmax), n=len)
    numpal= length(leg)
    pal= rainbow(numpal, start=0, end=4/6)[numpal:1] # revert colors: blue=low, red=high
    # Draw legend
    legend("right", bty="n", fill=c(pal, clr_nodata),legend=c(leg,"NA"), title=title)
    # Update limits to be the same as in the legend
    zmin= min(leg)
    zmax= max(leg)
    # More palette entries for data plotting
    numpal= 255
    pal= rainbow(numpal, start=0, end=4/6)[numpal:1] # revert colors: blue=low, red=high
    # Data
    for (ir in 1:nr) {
      clr= rep(clr_nodata, nc)
      inds= which(c(grid$matrix[ir,]) != grid$nodata_value)
      if (length(inds) > 0) {
        if (zmax == zmin) {
          clr[inds]= pal[1]
        } else {
          clr[inds]= pal[round((as.double(grid$matrix[ir,])[inds]-zmin)/(zmax-zmin)*(numpal-1)+1)]
        }
      }
      rect(xleft=cell_xmin(1:nc,xll,cz), xright=cell_xmax(1:nc,xll,cz),
        ybottom=cell_ymin(ir,nr,yll,cz), ytop=cell_ymax(ir,nr,yll,cz),
        border=NA, col=clr)
    }
    # Draw additional objects/annotations
    fun(...)
  } else {
    clr_nodata= "white"
    z= unique(as.numeric(grid$matrix))
    f= approxfun(x=z, y=1:length(z), method="constant", rule=1)
    pal= rainbow(length(z), start=0, end=4/6)[round(runif(length(z))*(length(z)-1)+1)]
    pal[z == grid$nodata_value]= clr_nodata
    # Data
    for (ir in 1:nr) {
      clr= pal[f(grid$matrix[ir,])]
      rect(xleft=cell_xmin(1:nc,xll,cz), xright=cell_xmax(1:nc,xll,cz),
        ybottom=cell_ymin(ir,nr,yll,cz), ytop=cell_ymax(ir,nr,yll,cz),
        border=NA, col=clr)
    }
    # Draw additional objects/annotations
    fun(...)
    # Legend
    legend("right", bty="n", fill=c("red", clr_nodata),legend=c("data","NA"))
  }
  # Draw outline
  rect(xleft=xll, ybottom=yll, xright=xll+nc*cz, ytop=yll+nr*cz, col=NA)
  # Save to file
  if (!is.null(file)) {
    if (file.exists(file) && (!replace)) stop(paste("Output file '",file,"' already exists.",sep=""))
    savePlot(filename=file, type=c(type), device = dev.cur())
    graphics.off()
  }
  # Reset and return
  par(op)
  return(invisible(NULL))
}

################################################################################

#' Lenght of a line in 2D
#'
#' Returns the length of a line, represented as a vector of x,y-coordinates. The
#' vector is expected to be of type character and each element is interpreted as
#' a pair of x and y coordinates. Merging the x and y info into a single value
#' allows for use of the function as the 'FUN' argument of \code{tapply}.
#'
#' @param xy_strings A vector of positions (as strings). Each position consists
#'   of the x and y coordinate separated by a single separator character.
#' @param sepchar The character separating the x and y value.
#'
#' @return The line's length.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @examples
#' # A data frame containing two line features
#' lineSet= data.frame(
#'   id= c(1,1,2,2,2),
#'   coords= c("0,0","1,1","2,2","1,2","1,1"),
#'   stringsAsFactors=FALSE
#' )
#' L= tapply(X=lineSet$coords, INDEX=lineSet$id, FUN=lineLength, sepchar=",")
#' print(L)

lineLength= function(xy_strings,sepchar) {
  tmp= unlist(strsplit(xy_strings, split=sepchar, fixed=TRUE))
  x= as.numeric(tmp[seq(from=1, to=length(tmp)-1, by=2)])
  y= as.numeric(tmp[seq(from=2, to=length(tmp), by=2)])
  return( sum((diff(x)^2 + diff(y)^2)^0.5) )
}

################################################################################

#' Zonal statistics for a grid of continuous data
#'
#' The function expects two input grids, one with data of a continuous variable
#' and another one with integer codes defining zones (classes). For each unique
#' zone, the (spatially) corresponding values of the continuous variable are
#' analyzed and a statistics table is returned. The function can be used, for
#' example, to compute an elevation statistics for river basins.
#'
#' @param grid_zones A geogrid object defining the zones. The values in this
#'   grid should be integers.
#' @param grid_data A geogrid object holding data of a continuous variable
#'   (like elevation, for example). The extend and resolution of this grid must
#'   be the same as for \code{grid_zones}.
#' @param minCoverage A numeric value specifying the minimum coverage of the
#'   valid grid cells in \code{grid_zones} by the valid values in
#'   \code{grid_data}. The value must be in range \eqn{0 < minCoverage <= 1}.
#'   If the actual coverage is less than the value of \code{minCoverage}, a
#'   warning is generated. Note that the coverage is analyzed globally 
#'   but not for the individual zones.
#' @param prefix A character string used as a prefix when creating column names
#'   for the output table. Must be compatible with R's convention for names.
#'   See the return value for details.
#' @param addMedian Should the median be added to the statistics? (logical)
#' @param addQuartiles Should quartiles be added to the statistics? (logical)
#' @param addExtremes Should min and max be added to the statistics? (logical)
#' @param sdigits The number of significant digits in the output table. Should
#'   be an integer >= 1.
#' @param silent Print diagnostic messages? (logical)
#'
#' @return A data frame with at least two columns 'id' and 'avg'. The 'id'
#'   columns holds the zones' IDs (unique values in \code{grid_zone}). The 'avg'
#'   column contains the arithmetic mean of the spatially corresponding values
#'   in \code{grid_data}. If \code{addMedian} is \code{TRUE}, a column 'med' is
#'   appended. If \code{addQuartiles} is \code{TRUE}, the additional columns 'q25'
#'   and 'q75' are present which hold the values of the lower and upper quartiles.
#'   If \code{addExtremes} is \code{TRUE}, there are two additional
#'   columns named 'min' and 'max'.
#'   If the value of \code{prefix} is a non-empty string, this string is used as
#'   a prefix to the above-mentioned column names.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

geogrid.zones.continuous = function(grid_zones, grid_data, minCoverage=1,
  prefix="", addMedian=FALSE, addQuartiles=FALSE, addExtremes=FALSE, sdigits=6, silent=TRUE
) {
  if (!is.geogrid(grid_zones)) stop("Argument 'grid_zones' is not a valid geogrid object.")
  if (!is.geogrid(grid_data)) stop("Argument 'grid_data' is not a valid geogrid object.")
  checkArg(arg=minCoverage, len=1, type="numeric")
  if ((minCoverage <= 0) || (minCoverage > 1))
    stop("Argument 'minCoverage' must be > 0 and <= 1.")
  checkArg(arg=prefix, len=1, type="character")
  checkArg(arg=addMedian, len=1, type="logical")
  checkArg(arg=addQuartiles, len=1, type="logical")
  checkArg(arg=addExtremes, len=1, type="logical")
  checkArg(arg=sdigits, len=1, type="integer")
  if (sdigits < 0) sdigits= 1
  checkArg(arg=silent, len=1, type="logical")
  # Check compatibility of input grids
  if (!geogrid.compatible(grid_zones,grid_data))
    stop("Input grids 'grid_zones' and 'grid_data' are not spatially compatible.")
  # Check coverage
  if (!silent) print("Analyzing data coverage of zones...")
  n= sum(grid_zones$matrix != grid_zones$nodata_value)
  if (n == 0)
    stop("The grid defining the zones does not contain any valid data.")
  notCovered= which((grid_data$matrix == grid_data$nodata_value) &
    (grid_zones$matrix != grid_zones$nodata_value))
  coverage= (1 - length(notCovered) / n)
  if (coverage < minCoverage)
    warning(paste("Coverage of zones by data is ",signif(coverage,5),"% but",
      " minimum coverage of ",minCoverage," has been requested.",sep=""))
  if (coverage == 0)
    stop("Coverage of zones by data is zero.")
  # Fill gaps if some exist
  if (length(notCovered) > 0) {
    marker= min(c(grid_data$matrix,grid_data$nodata_value)) - 1  # get unused marker value
    grid_data$matrix[notCovered]= marker                         # mark cells to fill
    if (!silent) print("Filling uncovered cells...")
    tmp= geogrid.fillGaps(grid=grid_data, removevalue=marker, maxdist=NULL,
      removeAll=TRUE, silent=silent)
    grid_data= tmp
    rm(tmp,marker)
  }
  # Allocate result table
  info= data.frame(id= unique(as.numeric(grid_zones$matrix[grid_zones$matrix !=
    grid_zones$nodata_value])))
  # Compute requested statistic for zones
  # (1) Average : Always computed
  if (!silent) print("Computing zonal averages...")
  tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=mean)
  tmp= data.frame(id=names(tmp), new=signif(as.numeric(tmp),sdigits))
  names(tmp)[names(tmp) == "new"]= paste(prefix,"avg",sep="")
  tmp= subset(tmp, tmp$id != grid_zones$nodata_value)
  info= merge(x= info, y=tmp, by="id", all=TRUE)
  rm(tmp)
  if (any(is.na(info[,])))
    stop("Failed to compute statistics for zones.")
  # (2) Median
  if (addMedian) {
    if (!silent) print("Computing zonal medians...")
    tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=median)
    tmp= data.frame(id=names(tmp), new=signif(as.numeric(tmp),sdigits))
    names(tmp)[names(tmp) == "new"]= paste(prefix,"med",sep="")
    tmp= subset(tmp, tmp$id != grid_zones$nodata_value)
    info= merge(x= info, y=tmp, by="id", all=TRUE)
    rm(tmp)
  }
  # (3) Quartiles
  if (addQuartiles) {
    if (!silent) print("Computing zonal quartiles...")
    for (p in c(0.25,0.75)) {
      tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=quantile, probs=c(p))
      tmp= data.frame(id=names(tmp), new=signif(as.numeric(tmp),sdigits))
      names(tmp)[names(tmp) == "new"]= paste(prefix,"q",round(p*100),sep="")
      tmp= subset(tmp, tmp$id != grid_zones$nodata_value)
      info= merge(x= info, y=tmp, by="id", all=TRUE)
      rm(tmp)
      if (any(is.na(info[,])))
        stop("Failed to compute statistics for zones.")
    }
  }
  # (4) Extremes
  if (addExtremes) {
    if (!silent) print("Computing zonal extremes...")
    for (fun in c("min","max")) {
      if (fun=="min") tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=min)
      if (fun=="max") tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=max)
      tmp= data.frame(id=names(tmp), new=signif(as.numeric(tmp),sdigits))
      names(tmp)[names(tmp) == "new"]= paste(prefix,fun,sep="")
      tmp= subset(tmp, tmp$id != grid_zones$nodata_value)
      info= merge(x= info, y=tmp, by="id", all=TRUE)
      rm(tmp)
      if (any(is.na(info[,])))
        stop("Failed to compute statistics for zones.")
    }
  }
  # Check names in table
  if (!identical(names(info),make.names(names(info))))
    stop("Bad column names in result table. Check the 'prefix' argument.")
  # Return table
  return(info)
}

################################################################################

#' Zonal statistics for a grid of classified data
#'
#' The function expects two input grids, one with integers representing classes
#' and another one with integer codes defining zones (classes). For each unique
#' zone, the shares of the (spatially) corresponding classes computed. The
#' function can be used, for example, to compute the areal fractions of
#' different land use classes for river basins.
#'
#' @param grid_zones A geogrid object defining the zones. The values in this
#'   grid should be integers.
#' @param grid_data A geogrid object whose data represent classes
#'   (of land use, for example). The extend and resolution of this grid must
#'   be the same as for \code{grid_zones}.
#' @param minCoverage A numeric value specifying the minimum coverage of the
#'   valid grid cells in \code{grid_zones} by the valid values in
#'   \code{grid_data}. The value must be in range \eqn{0 < minCoverage <= 1}.
#'   If the actual coverage is less than the value of \code{minCoverage}, a
#'   warning is generated. Note that the coverage is analyzed globally 
#'   but not for the individual zones.
#' @param minShare A threshold value >= 0 and <= 1. If the
#'   areal share of a class in a particular zone is less than this, the areal
#'   share is set to zero, i.e. the class is 'dissolved'. However, this is not
#'   done if (1) the class is listed in the argument \code{keepClasses} or (2)
#'   the class is the one with the largest areal share in the particular zone.
#'   One can set \code{minShare} to 1 in order to keep the dominating class only.
#' @param keepClasses A vector of integers representing classes in the input
#'   grid \code{grid_data}. The areal share of those classes is never set to
#'   zero even if it is less than \code{minShare}.
#' @param prefix A character string used as a prefix when creating column names
#'   for the output table. Must not be empty and must be compatible with R's
#'   convention for names. See the return value for details.
#' @param reportArea Defaults to \code{TRUE}. If \code{FALSE}, the total area
#'   of the zones is omitted in the result table.
#' @param hideZeroFields Defaults to \code{FALSE}. If \code{TRUE}, columns
#'   containing nothing but zero-values are deleted from the result table.
#' @param ndigits The number of digits when rounding the results. Should
#'   be an integer >= 1.
#' @param silent Print diagnostic messages? (logical)
#'
#' @return A data frame with columns 'id' and as many columns as there are
#'   classes in \code{grid_data} unless \code{hideZeroFields} is \code{TRUE}.
#'   The latter colums hold the areal shares of
#'   the classes within the zones and the column names are created by
#'   concatenation of \code{prefix} and the class' integer code.
#'   If \code{reportArea} is \code{TRUE}, there will be an additional field
#'   'area' holding the total area covered by the zones. It can be used, for
#'   example, to re-convert the areal shares into areas.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

geogrid.zones.classified = function(grid_zones, grid_data, minCoverage=1,
  minShare=0, keepClasses=numeric(0), prefix="class_", reportArea=TRUE, hideZeroFields=FALSE, ndigits=2, silent=TRUE
) {
  if (!is.geogrid(grid_zones)) stop("Argument 'grid_zones' is not a valid geogrid object.")
  if (!is.geogrid(grid_data)) stop("Argument 'grid_data' is not a valid geogrid object.")

  checkArg(arg=minCoverage, len=1, type="numeric")
  if ((minCoverage <= 0) || (minCoverage > 1))
    stop("Argument 'minCoverage' must be > 0 and <= 1.")
  checkArg(arg=minShare, len=1, type="numeric")
  if ((minShare < 0) || (minShare > 1))
    stop("Argument 'minShare' must be >= 0 and =< 1.")
  checkArg(arg=keepClasses, len=NULL, type="numeric")
  checkArg(arg=prefix, len=1, type="character")
  if (nchar(prefix) == 0)
    stop("Argument 'prefix' must be a non-empty string.")
  checkArg(arg=reportArea, len=1, type="logical")
  checkArg(arg=ndigits, len=1, type="integer")
  if (ndigits < 0) ndigits= 2
  checkArg(arg=silent, len=1, type="logical")
  # Check compatibility of input grids
  if (!geogrid.compatible(grid_zones,grid_data))
    stop("Input grids 'grid_zones' and 'grid_data' are not spatially compatible.")
  # Check coverage
  if (!silent) print("Analyzing data coverage of zones...")
  n= sum(grid_zones$matrix != grid_zones$nodata_value)
  if (n == 0)
    stop("The grid defining the zones does not contain any valid data.")
  notCovered= which((grid_data$matrix == grid_data$nodata_value) &
    (grid_zones$matrix != grid_zones$nodata_value))
  coverage= (1 - length(notCovered) / n)
  if (coverage < minCoverage)
    warning(paste("Coverage of zones by data is ",signif(coverage,5),"% but",
      " minimum coverage of ",minCoverage," has been requested.",sep=""))
  if (coverage == 0)
    stop("Coverage of zones by data is zero.")
  # Fill gaps if some exist
  if (length(notCovered) > 0) {
    marker= min(c(grid_data$matrix,grid_data$nodata_value)) - 1  # get unused marker value
    grid_data$matrix[notCovered]= marker                         # mark cells to fill
    if (!silent) print("Filling uncovered cells...")
    tmp= geogrid.fillGaps(grid=grid_data, removevalue=marker, maxdist=NULL,
      removeAll=TRUE, silent=silent)
    grid_data= tmp
    rm(tmp,marker)
  }
  # Compute areas of the zones and use this as a basis for the result table
  if (!silent) print("Computing areas of zones...")
  info= tapply(X=grid_zones$matrix, INDEX=grid_zones$matrix, FUN=length) * (grid_zones$cellsize^2)
  info= data.frame(id=names(info), area=as.numeric(info))
  info= subset(info, info$id != grid_zones$nodata_value)
  # Compute areal shares
  classes= unique(grid_data$matrix[grid_data$matrix != grid_data$nodata_value])
  for (i in 1:length(classes)) {
    if (!silent) print(paste("Computing areal shares for class ",classes[i],"...",sep=""))
    colname= paste(prefix,classes[i],sep="")
    tmp= tapply(X=grid_data$matrix, INDEX=grid_zones$matrix, FUN=function(x,val) {sum(x == val)},
      val=classes[i]) * (grid_data$cellsize^2)
    tmp= data.frame(id=names(tmp), count=tmp)
    tmp= subset(tmp, tmp$id != grid_zones$nodata_value)
    names(tmp)[2]= colname
    info= merge(x=info, y=tmp, by="id", all=TRUE)
    if (any(is.na(info[,])))
      stop(paste("Failed to compute areal shares for class ",classes[i],".",sep=""))
    # Convert area to areal share
    info[,colname]= info[,colname]/info$area
  }
  # Optionally dissolve classes (i.e. set small shares to zero) but retain classes
  # in the specified set AND never dissolve the class with the largest share. The
  # latter condition guarantees that at least 1 class remains in each zone.
  if (!silent) print("Dissolving small areal fractions... ")
  if ((length(keepClasses) > 0) && (!all(keepClasses %in% classes)))
    stop("Vector 'keepClasses' contains class codes not being present in 'grid_data'.")
  largestShares= apply(X=info[,paste(prefix,classes,sep="")], MARGIN=1, FUN=max)
  for (i in 1:length(classes)) {
    colname= paste(prefix,classes[i],sep="")
    if (!(classes[i] %in% keepClasses)) {
      # Info for understanding: minShare is a scalar but largestShares is a vector
      info[(info[,colname] < minShare) & (info[,colname] < largestShares),colname]= 0
    }
  }
  # Rescale the areal shares of the remained classes to account for dissolved ones
  info[,paste(prefix,classes,sep="")]= info[,paste(prefix,classes,sep="")] *
    (1. / apply(X=info[,paste(prefix,classes,sep="")], MARGIN=1, FUN=sum))
  # Round to requested number of digits
  info[,paste(prefix,classes,sep="")]= round(info[,paste(prefix,classes,sep="")], ndigits)
  # Delete zero-only fields if requested
  if (hideZeroFields) {
	  for (i in 1:length(classes)) {
		colname= paste(prefix,classes[i],sep="")
		if (sum(info[,colname]) == 0) info[,colname]=NULL
	  }
  }
  # Check names in table
  if (!identical(names(info),make.names(names(info))))
    stop("Bad column names in result table. Check the 'prefix' argument.")
  # Keep area field?
  if (!reportArea) info$area=NULL
  # Return table
  return(info)
}

