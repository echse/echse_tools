
#' Cross-secions from digital elevation models
#'
#' The function extracts cross-sections from a digital elevation model. The
#' cross-sections are defined by their cutlines on the \emph{horizontal} earth's
#' surface. Multiple DEMs can be processed in a single call.
#'
#' @param filesDEM Vector of file names. Each element should point to a file
#'   containing elevation data in ESRI's ASCII grid format.
#' @param filesSHP Vector of file names. Each element should point to an ESRI
#'   shape file defining cutline(s) of cross-section(s).
#' @param minRelResol Factor to control the minimum horizontal resolution
#'   (offset spacing) for the extracted cross-sections. A value of 1.0 (default)
#'   sets the distance between adjacent offsets to the cell size of the DEM.
#'   Values of 0.5 or 2, for example, reduce or increase the spacing by a factor
#'   of 2, respectively.
#' @param id_field The ID field in the shape file's attribute table. The info
#'   in this fielsd is used to construct names for the extracted cross-sections.
#' @param zeroBase If \code{TRUE}, the elevation data for the individual
#'   cross-sections are transformed so that the elevation of the lowest point
#'   becomes zero. This facilitates the visual comparison of cross-sections. If
#'   \code{FALSE}, the elevation data as extracted from the DEM are retained.
#' @param ndigits Desired number of digits in \emph{horizontal} coordinates, i.e.
#'   offsets, x- and y-coordinates. The precision of vertical coordinates is
#'   determined by the original DEM data and is \emph{not} affected by this argument.
#' @param outdir Name of the directory for output files.
#' @param prefix Prefix used in the construction of file names. For each
#'   cross-section, the file name is obtained by concatenating this prefix with
#'   the cutline's ID. A suitable file name suffix is appended.
#' @param png Switch to turn the output of graphics files on (\code{TRUE}) or off.
#' @param replace If \code{TRUE}, existing files will be silently replaced.
#' @param silent If \code{TRUE}, some diagnostic messages are printed.
#'
#' @return The total number of extracted cross-sections (integer).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @note If a cutline consists of more than two points, the interior points are
#'   ignored and the straight connection of the end points is used as the cutline.
#'
#' @seealso \code{\link{geogrid.valuesAtPoints}}, \code{\link{lineShapeAsTable}}
#'
#' @export

xs.extractDEM= function(
  filesDEM,
  filesSHP,
  minRelResol=1.0,
  id_field= "id",
  zeroBase= TRUE,
  ndigits= 3,
  outdir= getwd(),
  prefix= "xs",
  png= FALSE,
  replace= FALSE,
  silent= TRUE
) {
  # Check args
  checkArg(arg=filesDEM, len=NULL, type="character")
  if (length(filesDEM) < 1)
    stop("Argument 'filesDEM' is empty.")
  checkArg(arg=filesSHP, len=NULL, type="character")
  if (length(filesSHP) < 1)
    stop("Argument 'filesSHP' is empty.")
  if (length(filesSHP) != length(filesDEM))
    stop("Arguments 'filesSHP' and 'filesDEM' differ in length.")
  checkArg(arg=minRelResol, len=1, type="numeric")
  if (minRelResol <= 0)
    stop("Argument 'minRelResol' must be > 0.")
  checkArg(arg=id_field, len=1, type="character")
  checkArg(arg=zeroBase, len=1, type="logical")
  checkArg(arg=ndigits, len=1, type="integer")
  ndigits= abs(ndigits)
  checkDir(outdir, mustBeEmpty=FALSE)
  checkArg(arg=prefix, len=1, type="character")
  checkArg(arg=png, len=1, type="logical")
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")
  # Compute
  for (i in 1:length(filesDEM)) {

    # Check files
    checkFileIn(filesDEM[i])
    checkFileIn(filesSHP[i])

    # Read files and perform basic checks
    if (!silent) print(paste("Reading elevation model '",basename(filesDEM[i]),"'...",sep=""))
    dem= geogrid.readAscii(filesDEM[i])
    if (!silent) print(paste("Reading shape file '",filesSHP[i],"'...",sep=""))
    cut= lineShapeAsTable(filesSHP[i], id_field=id_field, attribs=c(), endsOnly=TRUE)
    cut= cut$shp
    if (!(id_field %in% names(cut)))
      stop(paste("No field with name '",id_field,"' in shape file '",filesSHP[i],"'.",sep=""))
    # Create sampling points at cutlines
    if (!silent) print(paste("Creating sampling points from '",basename(filesSHP[i]),"'...",sep=""))
    pts=data.frame(id=c(), x=c(), y=c(), offset=c())
    minResolution= dem$cellsize * minRelResol
    for (k in 1:nrow(cut)) {
      len= sqrt((cut$x1[k]-cut$x2[k])^2 + (cut$y1[k]-cut$y2[k])^2)
      npts= ceiling(len/minResolution) + 1
      x= seq(from=cut$x1[k], to=cut$x2[k], length.out=npts)
      y= seq(from=cut$y1[k], to=cut$y2[k], length.out=npts)
      offset= seq(from=0, to=len, length.out=npts)
      pts= rbind(pts, data.frame(id=rep(cut[k,id_field],npts),
        x=round(x,ndigits), y=round(y,ndigits), offset=round(offset,ndigits)))
    }

    # Exctract data
    if (!silent) print(paste("Extracting data from DEM '",filesDEM[i],"'...",sep=""))
    pts$elevation= geogrid.valuesAtPoints(grid=dem, xy_table=pts, silent=silent)

    # Collect results
    if (!silent) print("Collecting results...")
    if (i == 1) {
      sav= pts
    } else {
      dup= which(unique(pts$id) %in% sav$id)
      if (any(dup))
        stop(paste("Cut line IDs are not unique in the input shape files.",
          " Duplicates are '",paste(unique(pts$id)[dup],collapse="', '"),"'.",sep=""))
      sav= rbind(sav, pts)
    }
  } # End of loop over input files

  # Write to files
  if (!silent) print("Creating output files...")
  for (id in unique(sav$id)) {

    # Select data for particular x-section
    thisXS=sav[sav$id==id,]

    # Transform elevations if requested
    if (zeroBase)
      thisXS$elevation= thisXS$elevation - min(thisXS$elevation)

    # Create text file
    file_out=paste(outdir,"/",prefix,format(id,scientific=FALSE),".txt",sep="")
    if (file.exists(file_out) && (!replace))
      stop(paste("File '",file_out,"' already exists.",sep=""))
    write.table(thisXS,file=file_out,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

    # Create graphics
    if (png) {
      file_out=paste(outdir,"/",prefix,format(id,scientific=FALSE),".png",sep="")
      if (file.exists(file_out) && (!replace))
        stop(paste("File '",file_out,"' already exists.",sep=""))
      plot(thisXS$offset, thisXS$elevation, type="l")
      legend("topright", bty="n", paste("ID:",format(id,scientific=FALSE)))
      savePlot(filename=file_out, type="png", device=dev.cur())
    }
  }
  return(length(unique(sav$id)))
}


