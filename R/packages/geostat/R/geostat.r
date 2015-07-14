#' Utilities for spatial interpolation
#'
#' Type \code{help(package="geostat")} to inspect the package description.
#'
#' @name geostat-package
#' @aliases geostat
#' @docType package
{}


# Wrappers for C/C++ library functions follow:

################################################################################

#' Distance of two points on a plane
#' 
#' Function to compute the distance of two points on a plane
#'
#' @param x0 x-coordinate of 1st point
#' @param y0 y-coordinate of 1st point
#' @param x1 x-coordinate of 2nd point
#' @param y1 y-coordinate of 2nd point
#'
#' @return The distance between the two points.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib geostat
#'
#' @examples
#' dist2d(0,0,1,1)   # sqrt(2)

dist2d= function(x0, y0, x1, y1) {
  if (any(c(length(x0),length(y0),length(x1),length(y0)) != 1)) {
    stop("Coordinates must be scalars.")
  }
  if (any(!is.numeric(c(x0,y0,x1,y0)))) {
    stop("Coordinates must be numeric.")
  }
  result= .C("gs_dist", as.double(x0), as.double(y0), as.double(x1),
    as.double(y1), res=double(1))$res
  return(result)
}

################################################################################

#' Point of the compass when looking from one point to another point
#' 
#' When looking from a point 'A' to another point 'B', the function computes the
#' point of the compass corresponding to 'B'. The result is the angle against
#' North (0 == North, 90 == East, 180 == South, 270 == West).
#'
#' @param x_from x-coordinate of the base point
#' @param y_from y-coordinate of the base point
#' @param x_to x-coordinate of the point located from the base point
#' @param y_to y-coordinate of the point located from the base point
#'
#' @return An angle in range 0...360.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib geostat
#'
#' @examples
#' northangle(0,0,1,1)    # 45
#' northangle(0,0,-1,-1)  # 225

northangle= function(x_from, y_from, x_to, y_to) {
  if (any(c(length(x_from),length(y_from),length(x_to),length(y_to)) != 1)) {
    stop("Coordinates must be scalars.")
  }
  if (any(!is.numeric(c(x_from,y_from,x_to,y_to)))) {
    stop("Coordinates must be numeric.")
  }
  result= .C("gs_northangle", as.double(x_from), as.double(y_from),
    as.double(x_to), as.double(y_to), res=double(1))$res
  return(result)
}

################################################################################

#' Sector where a point is located when looking from another point
#' 
#' When looking from a point 'A' to another point 'B', the function computes the
#' sector in which point 'B' is located. A sector is simply a range of angles.
#' Sector indices increase in clock-wise direction.
#' For example, if the number of sectors in 4, and the origin is at 0 degrees,
#' sector 1 would cover the range from 0 to 90 degrees (North-East sector).
#'
#' @param nsectors The number of sectors to use (minimum 1). Each sector will
#'   be 360 / \code{nsectors} degrees wide. Setting \code{nsectors} to 1 results
#'   in a nearest-neighbor approach (also known as Thiessen method).
#' @param angle_of_origin Defines the angle of origin, i.e. the boundary between
#'   the first and the last sector. If \code{angle_of_origin} is zero, the
#'   boundary between sector 1 and sector \code{nsectors} is at zero degrees
#'   (i.e. the North-South axis). If, for example, \code{nsectors} is 4 and 
#'   \code{angle_of_origin} is 90, the sector with index 1 covers the range from
#'   90 to 180 degrees.
#' @inheritParams northangle
#'
#' @return The index of the sector as an integer in range 1 to \code{nsectors}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib geostat
#'
#' @examples
#' nsources= 1000
#' xrng= c(0, 100)
#' yrng= c(0, 100)
#' x_tar= 50
#' y_tar= 50
#' nsectors= 4
#' orig_angle= 45.
#' src= data.frame(x= min(xrng) + runif(nsources) * diff(xrng),
#'   y= min(yrng) + runif(nsources) * diff(yrng))
#' src$sect= NA
#' for (i in 1:nrow(src))
#'   src$sect[i]= sector(nsectors,orig_angle,x_tar,y_tar,src$x[i],src$y[i])
#' plot(xrng, yrng, type="n", xlab="x", ylab="y", asp=1)
#' clr= rainbow(nsectors)
#' for (i in 1:nsectors)
#'   points(src$x[src$sect==i],src$y[src$sect==i],pch=i,cex=0.5,col=clr[i])

sector= function(nsectors, angle_of_origin, x_from, y_from, x_to, y_to) {
  if (any(c(length(x_from),length(y_from),length(x_to),length(y_to)) != 1)) {
    stop("Coordinates must be scalars.")
  }
  if (any(!is.numeric(c(x_from,y_from,x_to,y_to)))) {
    stop("Coordinates must be numeric.")
  }
  if ((!is.numeric(angle_of_origin)) || (length(angle_of_origin) != 1)) {
    stop("Angle of origin must be a scalar number.")
  }
  if ((!is.numeric(nsectors)) || (length(nsectors) != 1)) {
    stop("Number of sectors must be a scalar integer.")
  }
  result= .C("gs_sector", as.integer(nsectors), as.double(angle_of_origin),
    as.double(x_from), as.double(y_from), as.double(x_to), as.double(y_to),
    res=integer(1))$res
  return(result)
}

################################################################################

#' Find suitable neighboring points for spatial interpolation by sector search
#' 
#' Given a point 'T' with coordinates \code{x_tar, y_tar} and a set of points
#' 'S' defined by the coordinate vectors \code{x_src} and \code{y_src}, the
#' function identifies those points among the set 'S' which are suitable for
#' estimating the an unknown value of a variable at point 'T' from the values at
#' the points 'S'. This is basically the initial step of every spatial
#' interpolation procedure.
#' The function uses a sector search approach, i.e. the surrounding of point 'T'
#' is divided into a user-specified number of sectors. For each sector, a
#' \emph{single} point in the set 'S' is identified. This is the point in that
#' particular sector which is closest to point 'T'. A sector can be empty.
#'
#' @param x_tar x-coordinate of the target point (scalar)
#' @param y_tar y-coordinate of the target point (scalar)
#' @param x_src x-coordinates of the source points (vector)
#' @param y_src y-coordinates of the source points (vector)
#' @param norigins The number of sector origins to be testet. If
#'   \code{norigins} is > 1, the search sectors will be rotated with an
#'   increment of 90 / \code{norigins} degrees around the center point given by
#'   \code{x_tar} and \code{y_tar}. The returned results correspond to the
#'   best angle of origin tested. The criterias of optimality are the number of
#'   non-empty sectors and the average distance between the target point and the
#'   source points. If \code{nsectors} == 1, sector rotation does not make
#'   sense and \code{norigins} = 1 should be used. In a standard setting with
#'   \code{nsectors} ==4, one would typically set \code{norigins} to a value in
#'   range 1 through 4. Higher values may give better results at the cost of
#'   increased computation time.
#' @inheritParams sector
#'
#' @return A list with the following two components:
#'   \item{selected}{A logical vector or the same length as \code{x_src}. The
#'     \code{TRUE} elements identify the selected source points. The order of
#'     elements is the same as in \code{x_src} and \code{y_src}.}
#'   \item{distance}{A numeric vector or the same length as \code{x_src}
#'      containing the distances between the target point and all source points.
#'      The order of elements is the same as in \code{x_src} and \code{y_src}.}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib geostat
#'
#' @examples
#' ndata= 16
#' nfigcol= 4
#' xrng= c(0, 100)
#' yrng= c(0, 100)
#' nsectors= 8
#' norigins= 10
#' data= data.frame(
#'   x= min(xrng) + runif(ndata) * diff(xrng),
#'   y= min(yrng) + runif(ndata) * diff(yrng)
#' )
#' old.par= par(mar=c(1,1,1,1))
#' split.screen(c(ceiling(ndata/nfigcol),nfigcol))
#' for (i in 1:nrow(data)) {
#'   x_tar= data$x[i]
#'   y_tar= data$y[i]
#'   x_src= data$x[1:nrow(data) != i]
#'   y_src= data$y[1:nrow(data) != i]
#'   nb= neighbors(x_tar, y_tar, x_src, y_src, nsectors, norigins)
#'   inds= which(nb$selected)
#'   avgdist= mean(nb$dist[nb$selected])
#'   screen(i)
#'   plot(data$x,data$y,type="n",xlab="x",ylab="y", xaxt="n", yaxt="n" ,asp=1)
#'   points(x_tar,y_tar,pch=1)
#'   points(x_src,y_src,pch=3,col="grey")
#'   for (k in 1:length(inds)) {
#'     lines(c(x_tar,x_src[inds[k]]), c(y_tar,y_src[inds[k]]), col="blue")
#'   }
#'   mtext(side=3,paste(round(avgdist,3)," (",sum(nb$selected),")",sep=""))
#' }
#' close.screen(all=TRUE)
#' par(old.par)

neighbors= function(x_tar, y_tar, x_src, y_src, nsectors, norigins) {
  if (any(c(length(x_tar),length(y_tar)) != 1)) {
    stop("Target coordinates must be scalars.")
  }
  if (any(!is.numeric(c(x_tar,y_tar)))) {
    stop("Target coordinates must be numeric.")
  }
  if (any(c(length(x_src),length(y_src)) == 0)) {
    stop("Source coordinates must be vectors.")
  }
  if (any(!is.numeric(c(x_src,y_src)))) {
    stop("Source coordinates must be numeric.")
  }
  if (length(x_src) != length(y_src)) {
    stop("Source coordinates vectors differ in length.")
  }
  if ((!is.numeric(nsectors)) || (length(nsectors) != 1)) {
    stop("Number of sectors must be a scalar integer.")
  }
  if ((!is.numeric(norigins)) || (length(norigins) != 1)) {
    stop("Number of sector origins must be a scalar integer.")
  }
  result= .C("gs_neighbors", as.double(x_tar), as.double(y_tar),
    as.integer(length(x_src)), as.double(x_src), as.double(y_src),
    as.integer(nsectors), as.integer(norigins),
    mask= integer(length(x_src)), dist= double(length(x_src)),
    status= integer(1))
  if (result$status != 0) {
    stop("Failed to select neighbors. Check input arguments.")
  }
  return(list(
    selected=as.logical(result$mask),
    distance=result$dist
  ))
}

################################################################################

#' Compute weights for inverse-distance spatial interpolation
#' 
#' Given a point 'T' with coordinates \code{x_tar, y_tar} and a set of points
#' 'S' defined by the coordinate vectors \code{x_src} and \code{y_src}, the
#' function computes interpolation weights for all point in 'S'. These weights
#' can be used to estimate the an unknown value of a variable at point 'T' from
#' known values at the points 'S'. The weights are computed by the
#' inverse-distance approach, i.e. the value of the weight decreases as the
#' distance between 'T' and the particular point from set 'S' increases.
#'
#' @param power The power to be used when computing inverse distance weights.
#'   In many applications, a value of 2 is used as a default. Higher values lead
#'   to increased weights for nearer points. A value of zero makes the weights
#'   independend of the points' distance and results in simple averaging.
#' @inheritParams neighbors
#'
#' @return A numeric vector of the same length as \code{x_src} holding the
#'   inverse-distance interpolation weights. The order of elements is the same
#'   as in \code{x_src} and \code{y_src}.
#'
#' @note The selection of suitable points from set
#'   'S' is performed by the \code{\link{neighbors}} method. A weight of zero is
#'   assigned to all points from 'S' not being part of the selection (because
#'   preference was given to a nearer point in the same sector).
#'
#'   Use of the return value: The returned vector can be used to filter
#'   the set of source points. For that purpose, the return vector needs to be
#'   converted into a logical mask by checking which elements differ
#'   significantly from zero using something like \code{mask = x > 0.01}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#' @useDynLib geostat
#'
#' @examples
#' # Demonstrates the regionalization of point observations to a regular grid
#' nsources= 16
#' xrng= c(0,100)
#' yrng= c(0,100)
#' dx= 4
#' dy= 4
#' nsectors= 12
#' norigins= 10
#' power= 2.
#' ncolors= 100
#' # Generate some source data
#' sources= data.frame(
#'   x= min(xrng) + runif(nsources) * diff(xrng),
#'   y= min(yrng) + runif(nsources) * diff(yrng),
#'   value= rnorm(mean=0, sd=1, n=nsources)
#' )
#' # Generate target locations (regular grid)
#' nx= diff(xrng)/dx
#' ny= diff(yrng)/dy  
#' targets= data.frame(
#'   x=rep(seq(from=min(xrng)+dx/2, to=max(xrng)-dx/2, length.out=nx), ny),
#'   y=sort(rep(seq(from=min(yrng)+dy/2, to=max(yrng)-dy/2, length.out=ny), nx))
#' )
#' plot(xrng,yrng,type="n",xlab="",ylab="", xaxt="n", yaxt="n" ,asp=1, bty="n")
#' clr= terrain.colors(ncolors)
#' # Show interpolation results (target locations)
#' for (i in 1:nrow(targets)) {
#'   weights= idweights(targets$x[i], targets$y[i], sources$x, sources$y,
#'     nsectors, norigins, power)
#'   if (abs(sum(weights)-1) > 0.0001) stop("Sum of weights not 1.")
#'   val= sum(weights * sources$value)
#'   clrind= round((val-min(sources$value))/
#'   (max(sources$value)-min(sources$value))*(ncolors-1))+1
#'   rect(xleft=targets$x[i]-dx,xright=targets$x[i]+dx,
#'     ybottom=targets$y[i]-dx,ytop=targets$y[i]+dx,col=clr[clrind],border=NA)
#' }
#' # Show original values (source locations)
#' vect_clrinds= round((sources$value-min(sources$value))/
#'   (max(sources$value)-min(sources$value))*(ncolors-1))+1
#' points(sources$x,sources$y,pch=20,cex=2,col=clr[vect_clrinds])
#' points(sources$x,sources$y,pch=1,cex=2)

idweights= function(x_tar, y_tar, x_src, y_src, nsectors=4, norigins=1, power=2) {
  if (any(c(length(x_tar),length(y_tar)) != 1)) {
    stop("Target coordinates must be scalars.")
  }
  if (any(!is.numeric(c(x_tar,y_tar)))) {
    stop("Target coordinates must be numeric.")
  }
  if (any(c(length(x_src),length(y_src)) == 0)) {
    stop("Source coordinates must be vectors.")
  }
  if (any(!is.numeric(c(x_src,y_src)))) {
    stop("Source coordinates must be numeric.")
  }
  if (length(x_src) != length(y_src)) {
    stop("Source coordinates vectors differ in length.")
  }
  if ((!is.numeric(nsectors)) || (length(nsectors) != 1)) {
    stop("Number of sectors must be a scalar integer.")
  }
  if ((!is.numeric(norigins)) || (length(norigins) != 1)) {
    stop("Number of sector origins must be a scalar integer.")
  }
  if ((!is.numeric(power)) || (length(power) != 1)) {
    stop("Power must be a scalar numeric value.")
  }

  result= .C("gs_idweights", as.double(x_tar), as.double(y_tar),
    as.integer(length(x_src)), as.double(x_src), as.double(y_src),
    as.integer(nsectors), as.integer(norigins), as.double(power),
    weights= double(length(x_src)), status= integer(1))
  if (result$status != 0) {
    stop("Failed to compute weights. Check input arguments.")
  }
  return(result$weights)
}

################################################################################

#' Generates an external input locations table for an ECHSE-based model
#' 
#' This method creates the so-called external input locations table for an
#' ECHSE-based model. This table is used to assign a (set of) station(s) and
#' weight(s) to all simulated objects for a particular variable. In hydrological
#' models, the table usually assigns at least a (set of) raingage(s) to all
#' sub-basins of a catchment.
#'
#' @param file_targets File holding a table of target locations. The table must
#'   have (at least) three columns with names 'id' (default), 'x', and 'y'. The field
#'   separator can be specified with the \code{colsep} argument. The 'x' and
#'   'y' field should contain numeric coordinates while the id-field contains
#'   unique identifiers for all target locations. In (semi)-distributed
#'   hydrological modeling, this is typically a table of the coordinates of the 
#'   sub-basins' centers of gravity.
#' @param files_sources A \emph{named} vector of file names. Each of the files
#'   is expected to hold a table of source locations for a particular variable.
#'   The names of the variables are inferred from the names of the vectors'
#'   elements. In the context of hydrological modeling, these files typically
#'   list the coordinates of rain gages and sensors of other meteorological
#'   variables. All files must be in the same format as the file supplied as
#'   \code{file_targets}. The id-fields contain the station names/numbers.
#' @param idField_targets Name of the field in \code{file_targets} containing
#'   ID strings for the locations. Defaults to 'id'.
#' @param idField_sources Name of the field in \code{file_sources} containing
#'   ID strings for the locations. Defaults to 'id'.
#' @param colsep Field separator used in all input files and the output file.
#'   Defaults to TAB.
#' @param nsectors See corresponding argument of the \code{\link{idweights}}
#'   method. The default value of 4 is often used and results in a so-called
#'   'quadrant search'.
#' @param norigins See corresponding argument of the \code{\link{idweights}}
#'   method.
#' @param power See corresponding argument of the \code{\link{idweights}}
#'   method. Defaults to 2.
#' @param file_result Name of the result file created by this function.
#' @param ndigits The number of digits to be used when printing the computed
#'   weights to \code{file_result}. This value is passed as the second argument
#'   to the \code{\link{round}} method. Records for those source locations whose 
#'   \emph{rounded} weight is zero are dropped when creating \code{file_result}.
#' @param overwrite Is it OK to overwrite an existing result file? Defaults to
#'   FALSE.
#' @inheritParams idweights
#'
#' @return \code{NULL}
#'
#' @note For more information on format and use of the result file please have a
#'   look at the ECHSE documentation (reference given below). In particular, one
#'   should read the information on external input variables and the
#'   configuration file's keyword 'table_externalInput_locations'.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Demonstrates the assignment of hydro-meteorological sensors to sub-basins
#' # Dummy table of sub-basin coordinates
#' file_basins= tempfile()
#' tbl= data.frame(id= c("sb1","sb2"), x=runif(2), y=runif(2))
#' write.table(x= tbl, file=file_basins, sep="\t", row.names=FALSE, quote=FALSE)
#' # Dummy table of rain gages
#' file_gages= tempfile()
#' tbl= data.frame(id= c("g1","g2","g3"), x=runif(3), y=runif(3))
#' write.table(x= tbl, file=file_gages, sep="\t", row.names=FALSE, quote=FALSE)
#' # Dummy table of temperature sensors
#' file_temps= tempfile()
#' tbl= data.frame(id= c("t1","t2","t3","t4"), x=runif(4), y=runif(4))
#' write.table(x= tbl, file=file_temps, sep="\t", row.names=FALSE, quote=FALSE)
#' # Compute and show result table
#' file_result= tempfile()
#' externalInputLocationsTable(
#'   file_targets= file_basins,
#'   files_sources= c(precip= file_gages, temper= file_temps),
#'   colsep="\t", nsectors= 4, norigins= 2, power= 2,
#'   file_result=file_result)
#' print(read.table(file=file_result))
#' file.remove(file_basins, file_gages, file_temps, file_result)
#'
#' @references Kneis, D.: Eco-Hydrological Simulation Environment (ECHSE)
#'   - Documentation of the Generic Components.

externalInputLocationsTable= function(
  file_targets,
  files_sources,
  idField_targets="id",
  idField_sources="id",
  colsep="\t",
  nsectors= 4, norigins= 4, power= 2,
  file_result, ndigits= 3, overwrite= FALSE
) {

  # Basic input checks
  if (any(names(files_sources) == ""))
    stop("Missing element name(s) in vector 'files_sources'.")

  # Read target locations
  if (!file.exists(file_targets))
    stop(paste("Table of target locations '",file_targets,"' not found.",sep=""))
  tar= read.table(file=file_targets, header=TRUE, sep=colsep)
  if (!all(c(idField_targets,"x","y") %in% names(tar))) {
    stop(paste("Missing column(s) in table of target locations '",file_targets,
    "'. Expected: '",idField_targets,"', 'x', 'y'.",sep=""))
  }
  tar[,idField_targets]= as.character(tar[,idField_targets])
  if (length(unique(tar[,idField_targets])) != nrow(tar))
    stop(paste("Non-unique entries in field '",idField_targets,"' of '",file_targets,"'.",sep=""))
    
  # Initialize result file
  if (file.exists(file_result) && (!overwrite))
    stop(paste("Output file '",file_result,"' already exists.",sep=""))
  ovect= c("object","variable","location","weight")
  write(ovect, file=file_result, ncolumns=length(ovect), sep=colsep, append=FALSE)

  # Loop through variables
  for (i in 1:length(files_sources)) {
    
    # Read source locations for this variable
    if (!file.exists(files_sources[i]))
      stop(paste("Table of source locations for variable '",names(files_sources)[i],
        "' ('",files_sources[i],"') not found.",sep=""))
    src= read.table(file=files_sources[i], header=TRUE, sep=colsep)
    if (!all(c(idField_sources,"x","y") %in% names(src))) {
      stop(paste("Missing column(s) in table of source locations for variable '",
        names(files_sources)[i],"' ('",
        files_sources[i],"'). Expected: '",idField_sources,"', 'x', 'y'.",sep=""))
    }
    src[,idField_sources]= as.character(src[,idField_sources])
    if (length(unique(src[,idField_sources])) != nrow(src))
      stop(paste("Non-unique entries in field '",idField_sources,"' of table of source locations",
        " for variable '", names(files_sources)[i],"' ('",files_sources[i],"').",sep=""))

    # Loop through target locations
    for (n in 1:nrow(tar)) {
      # Compute weights
      weights= idweights(tar$x[n], tar$y[n], src$x, src$y, nsectors, norigins, power)
      # Drop weights than would become zero after rounding
      inds= which(round(weights,ndigits) >= 1/(10^ndigits))
      for (k in 1:length(inds)) {
        ovect= c(tar[n,idField_targets], names(files_sources)[i] ,src[inds[k],idField_sources],
          round(weights[inds[k]],ndigits))
        write(ovect, file=file_result, ncolumns=length(ovect), sep=colsep, append=TRUE)
      }
    }
  }
  return(invisible(NULL))
}

################################################################################

#' Computes time series of residuals assuming a linear model
#'
#' The function takes a multi-location time series of an observed variable
#' \eqn{y} as input. Typically, \eqn{y} is a meteorological variable.
#' In addition, the function requires information on a static predictor
#' variable \eqn{z}, whose value is known at the locations where \eqn{y} has
#' been measured.
#' Assuming the existance of a \emph{linear} correlation \eqn{y= a * z + b}, the
#' function converts the original data (\eqn{y}) into (additive) residuals
#' \eqn{r = y - (a * z + b)}. The coefficients \eqn{a} and \eqn{b} are estimated
#' either for every single time step or time-averaged data.
#' The function returns both a time series of the residuals \eqn{r} as well as a
#' time series of the linear model's coefficients \eqn{a} and \eqn{b}.
#'
#' @param file_obs Input file containing the time series of observations. Column names
#'   are expected in the first row. The column names are interpreted as location
#'   IDs except for the one specified by \code{obs_colTime} (header of the time
#'   column). The strings in the time column can be in any convenient format,
#'   for example '2000-01-01 12:30:00'. The time column is copied to the output
#'   files 'as is'.
#' @param file_loc Input file containing a locations' attribute. Must be a table with
#'   two named columns at least, one holding the locations' IDs and the other
#'   one containing the value of a numerical attribute (typically elevation).
#'   For each location in \code{file_obs}, there must be one matching record.
#' @param file_resid Output file holding the time series of residuals. The
#'   file's layout is identical with that of the input file \code{file_obs}.
#' @param file_coeff Output file holding the time series of coefficients of the
#'   linear model. The first column in this file is identical with the time
#'   column from \code{file_obs}. The columns 'intercept' and 'slope' hold the
#'   model coefficients and the column 'r2' list the values of R-squared.
#' @param colsep Field separator used in all input and output files.
#'   Defaults to TAB.
#' @param obs_colTime Name of the column in \code{file_obs} containing time info
#'   instead of numerical data. Defaults to 'datetime'.
#' @param loc_colId Name of the column in \code{file_loc} holding the locations'
#'   IDs. Defaults to 'id'.
#' @param loc_colPred Name of the column in \code{file_loc} holding the values
#'   of the predictor variable. Defaults to 'z'.
#' @param average Logical argument (default \code{FALSE}). If \code{FALSE},
#'   regression coefficients are computed for individual records, i.e. for every
#'   single time step. If \code{TRUE}, the data are averaged over time before
#'   the coefficients of the linear model are estimated. These parameters are
#'   then applied to \emph{all} time steps. 
#' @param r2min Minimum value of R-squared for the linear model relating the
#'   observed variable to the predictor variable. The linear model is only used
#'   if the value of R-squared for a particular time step is >= \code{r2min}. In
#'   case of a lower value of R-squared, slope and intercept of the linear model
#'   are set to zero and the residuals become identical with the original data.
#'   The default for \code{r2min} is 0.36 which is equivalent to a coefficient
#'   of correlation of 0.6.
#' @param nsignif The number of \emph{significant} digits to be used when
#'   printing the residuals to \code{file_resid}. Defaults to 3.
#' @param overwrite Defaults to FALSE. If TRUE, existing output files are
#'   replaced without a warning.
#'
#' @return \code{NULL}
#'
#' @note When using the output of this function for spatial residual
#'   interpolation, one needs to take care of possible undesired extrapolation
#'   effects. Problems are likely to occur if the interpolated variable is (1) 
#'   physically limited to a certain range and (2) the value of the predictor
#'   variable at a target location is far outside the predictor variable's range
#'   at the source locations. Example: Rainfall is interpolated for a
#'   low-elevation site using rain gage data from higher sites. If the
#'   correlation between elevation and rainfall is strong but negative, the
#'   estimated rainfall for the low-elevation site may be negative.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Test input data
#' print("__________ Obs. data __________")
#' obs= data.frame(
#'   datetime=ISOdatetime(2000,1,1,0,0,0) + seq(from=0, to=2, by=1) * 3600,
#'   stationA=c(1,2,8), stationB=c(2,3,6), stationC=c(6,2,1), stationD=c(4,3,3))
#' print(obs)
#' write.table(obs, file=fObs<-tempfile(), col.names=TRUE,row.names=FALSE,sep="\t")
#' print("__________ Locations __________")
#' loc= data.frame(
#'   id= c("stationA","stationB","stationC","stationD"),
#'   elevation= c(2, 400, 1800, 900))
#' print(loc)
#' write.table(loc, file=fLoc<-tempfile(), col.names=TRUE,row.names=FALSE,sep="\t")
#' # Convert
#' residualTimeSeries(
#'   file_obs= fObs, file_loc= fLoc,
#'   file_resid=fResid<-tempfile(), file_coeff=fCoeff<-tempfile(),
#'   colsep="\t",
#'   obs_colTime="datetime", loc_colId="id", loc_colPred="elevation",
#'   nsignif=4, overwrite=FALSE
#' )
#' # Result of conversion
#' print("__________ Residuals __________")
#' resid= read.table(file=fResid, header=TRUE, sep="\t")
#' print(resid)
#' print("__________ Model fit __________")
#' coeff= read.table(file=fCoeff, header=TRUE, sep="\t")
#' print(coeff)
#' # Back conversion of residuals
#' print("__________ Back conv. __________")
#' back= resid
#' for (i in 2:ncol(resid))
#'   back[,i]= resid[,i] +
#'     (coeff$slope * loc$elevation[match(names(obs)[i],loc$id)] + coeff$intercept)
#' print(back)


residualTimeSeries= function(
  file_obs,
  file_loc,
  file_resid,
  file_coeff,
  colsep="\t",
  obs_colTime="datetime",
  loc_colId="id",
  loc_colPred="z",
  average= FALSE,
  r2min= 0.36,
  nsignif=3,
  overwrite=FALSE
) {
  # Read/check observations
  if (!file.exists(file_obs))
    stop(paste("File with observed time series '",file_obs,"' not found.",sep=""))
  obs= read.table(file=file_obs, header=TRUE, sep=colsep)
  if (!(obs_colTime %in% names(obs)))
    stop(paste("Missing column '",obs_colTime,"' in file '",file_obs,"'.",sep=""))
  if (ncol(obs) < 2)
    stop(paste("Too few columns in file '",file_obs,"'.",sep=""))
  obs[,obs_colTime]= as.character(obs[,obs_colTime])  # We don't need the time values

  # Read/check static predictor variable for locations
  if (!file.exists(file_loc))
    stop(paste("File with location attributes '",file_loc,"' not found.",sep=""))
  loc= read.table(file=file_loc, header=TRUE, sep=colsep)
  if (!all(c(loc_colId,loc_colPred) %in% names(loc)))
    stop(paste("Missing column(s) in file with location attributes '",file_loc,
    "'. Expecting at least '",loc_colId,"' and '",loc_colPred,"'.",sep=""))
  loc[,loc_colId]= as.character(loc[,loc_colId])
  if (length(unique(loc[,loc_colId])) != nrow(loc))
    stop(paste("Non-unique entries in '",loc_colId,"' column of '",file_loc,"'.",sep=""))

  # Identify data columns in time series input file
  obs_colsLocs= which(names(obs) != obs_colTime)

  # Get value of predictor and original name for locations in input time series
  inds= match(names(obs)[obs_colsLocs], make.names(loc[,loc_colId]), nomatch=NA)
  if (any(is.na(inds)))
    stop(paste("Missing record(s) in file '",file_loc,
      "' corresponding to locations(s) in file '",file_obs,"'.",
      " Info is lacking for locations(s): '",
      paste(names(obs)[obs_colsLocs][is.na(inds)],collapse="', '"),"'.",sep=""))
   predictors= loc[inds,loc_colPred]
   locnames= loc[inds,loc_colId]
   rm(inds)

  # Make sure that the values of the predictor variable are not constant
  # as this would lead to a singular case in model fitting
  if (length(unique(predictors)) == 1)
    stop("Value of the predictor variable is identical for all locations.")

  # Set up data frame with linear model parameters for all time steps
  # We don't use lm() here because this seems to be slow
  lmCoeffs= function(predicted,predictor) {
    if (length(predicted) < 2) {
      stop("Cannot fit linear model. Input vectors too short.")
    } else {
      var_predictor= var(predictor);
      var_predicted= var(predicted);
      if (var_predictor == 0) { # All z-values equal --> Infinite slope
        stop("Cannot fit linear model. Variance of predictor variable is zero.")
      } else {
        if (var_predicted == 0) { # All v-values equal --> Zero slope, zero correlation
          return(c(intercept=0, slope=0, r2=0)) # OK, but strange input
        } else { # The normal, well-behaved case
          covar= cov(predictor,predicted)
          slope= covar / var_predictor
          intercept= mean(predicted) - slope * mean(predictor)
          return(c(intercept= intercept, slope= slope,
            r2= (covar^2)/(var_predictor*var_predicted)))
        }
      }
    }
  }

  if (!average) { # Regression for every single time step
    pars= apply(X=obs[,obs_colsLocs], MARGIN=1, FUN=lmCoeffs, predictor=predictors)
    pars= as.data.frame(t(pars))
    pars= cbind(time=obs[,obs_colTime], pars)
    names(pars)[1]= obs_colTime
  } else { # Regression based on average data
    avg= apply(X=obs[,obs_colsLocs], MARGIN=2, FUN=mean)
    tmp= lmCoeffs(avg,predictors)
    pars= data.frame(time=obs[,obs_colTime], intercept=tmp[["intercept"]],
      slope=tmp[["slope"]], r2=tmp[["r2"]])
    names(pars)[1]= obs_colTime
    rm(tmp)
  }

  # Set coefficients for records with too low r2
  pars$slope[which(pars$r2 < r2min)]= 0
  pars$intercept[which(pars$r2 < r2min)]= 0

  # Apply linear model to convert original data into residuals
  for (i in 1:length(obs_colsLocs)) {
    obs[,obs_colsLocs[i]]= obs[,obs_colsLocs[i]] -
      (predictors[i] * pars$slope + pars$intercept)
  }

  # Create output file 1: Time series of linear model parameters
  if (file.exists(file_coeff) && (!overwrite))
    stop(paste("Output file '",file_coeff,"' already exists.",sep=""))
  pars$intercept= signif(pars$intercept, 6)
  pars$slope= signif(pars$slope, 6)
  pars$r2= round(pars$r2, 2)
  write.table(x=pars, file=file_coeff, sep=colsep, col.names=TRUE,
    row.names=FALSE, quote=FALSE)

  # Create output file 2: Time series of linear model parameters
  if (file.exists(file_resid) && (!overwrite))
    stop(paste("Output file '",file_resid,"' already exists.",sep=""))
  ovect= c(obs_colTime, locnames)
  write(x=ovect, ncolumns=length(ovect), file=file_resid, sep=colsep)
  obs[,obs_colsLocs]= signif(obs[,obs_colsLocs], nsignif)
  write.table(x=obs, file=file_resid, sep=colsep, col.names=FALSE,
    row.names=FALSE, quote=FALSE, append=TRUE)

  # Begin (check section)
  # Note: For checks, the number of digits must be chosen large enough
  if (FALSE) { # FALSE: disabled, TRUE: enabled
    resid= read.table(file=file_resid, header=TRUE, sep=colsep)
    coeff= read.table(file=file_coeff, header=TRUE, sep=colsep)
    for (i in 1:length(obs_colsLocs)) {
      resid[,obs_colsLocs[i]]= resid[,obs_colsLocs[i]] +
        (predictors[i] * coeff$slope + coeff$intercept)
    }
    obs= read.table(file=file_obs, header=TRUE, sep=colsep)
    for (i in 1:length(obs_colsLocs)) {
      maxerr= max(abs(resid[,obs_colsLocs[i]] - obs[,obs_colsLocs[i]]))
      if (maxerr > 10^(-nsignif+1))
        stop(paste("Error too large."))
    }
  }
  # End (check section)

  return(invisible(NULL))
}

