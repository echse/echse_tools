################################################################################

# Basic code for extraction of data from shape file found here:
#
# http://r-sig-geo.2731867.n2.nabble.com/How-to-extract-coordinates-values-from-a-shapefile-td5159868.html


#' Tabular representation of a shape file with line features
#'
#' The function converts an ESRI shape file with line features into a table of
#' coordinates. Feature attributes can be extracted as well.
#'
#' @param file Name/path of the shape file (with or without extension .shp). The
#'   file must contain \emph{line} features (not points or polygons). The lines
#'   must \emph{not} be multi-part lines, i.e. each line must be a continuous
#'   sequence of coordinates without breaks (gaps).
#' @param id_field Name of the column in the shape file's attribute table
#'   containing the unique object IDs.
#' @param attribs Vector of (further) attribute columns names. By default, no
#'   further attributes are extracted.
#' @param endsOnly If \code{TRUE}, the result table contains only the
#'   coordinates of the lines' end points but no intermediate coordinates. If
#'   \code{FALSE} (default), the full array of coordinates is returned.
#'
#' @return A list of two data frames with names 'shp' and 'att'. The latter
#'   data frame represents the shape file's attribute table. It contains the
#'   ID field defined by \code{id_field} and the attribute fields listed in the
#'   \code{attribs} argument or \code{NULL} if \code{attribs} is an empty vector.
#'   The contents of the 'shp' data frame depends on the value of argument
#'   \code{endsOnly}. See notes for details.
#'
#' @note If \code{endsOnly} is \code{FALSE}, the 'shp' data frame has the 3
#'   columns. The first column is the ID field with the name determined by
#'   \code{id_field}. The other two columns are 'x', and 'y'. The table
#'   contains the complete list of the lines' coordinates.
#'   The number of records in the result table is
#'   \eqn{\sum_{i=1}^n{L_i}} where \eqn{n} is the number of lines and \eqn{L_i}
#'   is the number of coordinates of a particular line with index \eqn{i}.
#'
#' If \code{endsOnly} is \code{FALSE}, the 'shp' data frame has 5 columns.
#' The first column is the ID field with the name determined by \code{id_field}.
#' The 4 remaining columns are 'x1', 'y1', 'x2', and 'y2'. The table contains
#' only the coordinates of the lines' start and end points in the fields 'x1', 'y1' and
#' 'x2', 'y2', respectively. The result table has as many records as there are
#' line features in the shape file.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @seealso \code{\link{slot}}, \code{\link{readShapeLines}} from package \code{maptools}
#'
#' @export
#'

lineShapeAsTable = function(file, id_field, attribs=c(), endsOnly=FALSE) {
  # Check input
  checkFileIn(file)
  checkArg(arg=id_field, len=1, type="character")
  if (length(attribs) != 0)
    checkArg(arg=attribs, len=NULL, type="character")
  checkArg(arg=endsOnly, len=1, type="logical")
  # Read shapefile
  x= maptools::readShapeLines(file)
  # Check existence/consistency of id_field and further attributes
  if (!(id_field %in% names(x@data)))
    stop(paste("ID column with name '",id_field,"' do(es) not exist in attribute table.",sep=""))
  if (length(attribs) > 0) {
    if (!all(attribs %in% names(x@data))) {
      stop("Column name(s) in 'attribs' do(es) not exist in attribute table.")
    }
    if (id_field %in% attribs)
      stop(paste("Column '",id_field,"' is already used as ID field and cannot be an attribute field too.",sep=""))
  }
  # Check ID column
  if (length(unique(x@data[,id_field])) != nrow(x@data))
    stop(paste("Duplicates detected in the attribute table's ID field '",id_field,"'.",sep=""))
  # Init temporary file to store coordinates
  f= tempfile()
  write(x= c(id_field,"x","y"), file=f, ncolumns=3, sep="\t")
  # Extract coordinates and attributes
  allLines = x@lines
  for (i in 1:length(allLines)) {
    thisLine = allLines[[i]]@Lines
    id= x@data[i,id_field]
    nSegm= length(thisLine)
    if (nSegm > 1)
      stop(paste("Line with ID '",id,"' is a multi-part line with ",nSegm," segments.",sep=""))
    for (j in 1:nSegm) {
      xy= sp::coordinates(thisLine[[j]])
      nCoords= nrow(xy)
      write.table(x= cbind(rep(id,nCoords), xy),
        file=f, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE)
    }
  }
  shp= read.table(f, header=TRUE, sep="\t")
  file.remove(f)
  # Transform into table of end-points, if requested
  if (endsOnly) {
    featureLens= tapply(X=shp[,id_field], INDEX=shp[,id_field], FUN=length)
    inds1= (cumsum(featureLens)-featureLens)+1
    inds2= cumsum(featureLens)
    shp= data.frame(id=shp[inds1,id_field], x1= shp[inds1,"x"], y1= shp[inds1,"y"], x2= shp[inds2,"x"], y2= shp[inds2,"y"])
    names(shp)[names(shp) == "id"]= id_field
  }
  if (length(attribs) > 0) {
    att= x@data[,c(id_field,attribs)]
  } else {
    att= NULL
  }
  return(list(shp=shp, att=att))
}

## TEST CODE
if (FALSE) {
  library(topocatch)
  setwd("/home/dkneis/progress/echse_proj/indMod/mahanadi/data/topocatch/out10")
  d = lineShapeAsTable("net.shp",id_field="id", attribs=c("class"), endsOnly=FALSE)
  plot(range(d$shp$x), range(d$shp$y), type="n")
  d$shp$xy= paste(d$shp$x,d$shp$y)
  lin= function(z) {
    zz=unlist(strsplit(z," ",fixed=TRUE))
    n= length(zz)
    lines(zz[seq(1,n-1,2)],zz[seq(2,n,2)])
    return(invisible(NULL))
  }
  nothing= tapply(X=d$shp$xy, INDEX=d$shp$id, FUN=lin)

  d2 = lineShapeAsTable("net.shp",id_field="id", attribs=c("class"), endsOnly=TRUE)
  for (i in 1:nrow(d2$shp)) {
    lines(x=c(d2$shp$x1[i],d2$shp$x2[i]), y=c(d2$shp$y1[i],d2$shp$y2[i]), col="red")
  }
}


