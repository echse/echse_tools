
# The workhorse function for generating a catchment grid.
# It is called by the pulic high-level function 'hydroData'.

catchments = function(
  fileDIR,                             # Input grid with flow direction codes
  fileSHP,                             # Input shape file with drainage lines
  fileCAT,                             # Output grid file with catchment IDs
  id_field="id",                       # ID field in the input shape file
  class_field="class",                 # Class field in the input shape file
  classes_with_catchment= c("rch"),   # Names of classes to which a catchment is to be assigned
  nbuffer=1,                           # Number of buffers in vector-raster conversion
  replace=FALSE,                       # Handling of existing output files
  silent=TRUE                          # Switch for diagnostic messages
) {
  # Check args
  # Bool
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")
  # Files
  checkFileIn(fileDIR)
  checkFileIn(fileSHP)
  checkFileOut(fileCAT, replace=replace)
  # Names
  checkArg(arg=id_field, len=1, type="character")
  checkArg(arg=class_field, len=1, type="character")
  checkArg(arg=classes_with_catchment, len=NULL, type="character")
  if (length(classes_with_catchment) < 1)
    stop("Empty vector of class names.")
  # Numeric
  checkArg(arg=nbuffer, len=1, type="integer")
  # Convert and filter shape file
  if (!silent) print("Converting shape file...")
  shptab= lineShapeAsTable(file=fileSHP, id_field=id_field,
    attribs=c(class_field), endsOnly=FALSE)
  shptab= merge(x=shptab$shp, y=shptab$att, by=id_field, all=TRUE)
  if (any(is.na(shptab)))
    stop("Failed to merge feature geometry and attributes.")
  if (!silent) print("Filtering features in shape file...")
  shptab= shptab[(shptab[,class_field] %in% classes_with_catchment),]
  if (nrow(shptab) == 0)
    stop(paste("Shape file '",fileSHP,"' does not contain objects with the",
      " specified class name(s) in field '",class_field,"' of the attribute",
      " table.",sep=""))
  # Read flow direction grid
  if (!silent) print("Reading flow direction grid...")
  fdir= geogrid.readAscii(fileDIR)
  # Set a save special code for cells touched by multiple lines (at junctions,
  # for example). This code must not be a valid line ID. The respective cells
  # are treated specially when building the catchments.
  codeConflict= min(c(shptab[,id_field])) - 1
  # Set a nodata value for the result grid
  nodata= min(codeConflict, -9998) - 1
  # Convert lines to grid (for use as a start grid in catchment building)
  if (!silent) print("Doing vector-to-raster conversion...")
  ini= linesToGrid(ncol=fdir$ncols, nrow=fdir$nrows, xllcorner=fdir$xllcorner,
    yllcorner=fdir$yllcorner, cellsize=fdir$cellsize, lineTable=shptab,
    code_undefined=nodata, code_conflict=codeConflict,
    nbuffer=nbuffer, silent=silent)
  # Build raw catchments
  if (!silent) print("Building raw catchments...")
  codeUndefined= min(codeConflict, nodata) - 1
  raw= buildCatchments(grid_init=ini, code_conflict=codeConflict,
    grid_flowdir=fdir, code_undefined=codeUndefined)
  # Exclude cells outside the extent of the flowdir grid
  raw$matrix[fdir$matrix == fdir$nodata_value]= raw$nodata_value
  # Fill undefined cells in the raw catchment grid
  if (!silent) print("Filling undefined cells...")
  cat= geogrid.fillGaps(grid=raw, removevalue=codeUndefined, maxdist=NULL,
    removeAll=TRUE, silent=silent)
  # Write result grid
  if (!silent) print("Writing result grid...")
  geogrid.writeAscii(grid=cat, file=fileCAT, replace=replace)
  return(invisible(NULL))
}

