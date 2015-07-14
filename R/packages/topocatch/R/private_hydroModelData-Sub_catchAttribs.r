# This function computes the basic catchment attributes (position and size).
# It is called by the pulic high-level function 'hydroData'.

catchAttribs= function(
  fileCAT,                   # Grid file with catchment codes
  prefix,                    # prefix for generation of catchment IDs frmo numerical codes
  silent                     # Switch for diagnostic messages
){
  # Check args
  checkFileIn(fileCAT)
  checkArg(arg=prefix, len=1, type="character")
  checkArg(arg=silent, len=1, type="logical")
  # Read grid
  if (!silent) print("Reading catchment grid...")  
  cat= geogrid.readAscii(fileCAT)
  # Init table
  info= data.frame(object= unique(as.numeric(cat$matrix[cat$matrix != cat$nodata_value])))
  # Compute positions (centers of gravity)
  if (!silent) print("Computing x-coordinates...")
  x_4row= cell_xavg(1:cat$ncols,cat$xllcorner,cat$cellsize)
  xcoords= matrix(x_4row, ncol=cat$ncols, nrow=cat$nrows, byrow=TRUE)
  tmp= tapply(X=xcoords, INDEX=cat$matrix, FUN=mean)
  tmp= data.frame(id=names(tmp), x=round(as.numeric(tmp)))
  tmp= subset(tmp, tmp$id != cat$nodata_value)
  info= merge(x= info, y=tmp, by.x="object", by.y="id", all=TRUE)
  rm(x_4row, xcoords, tmp)
  if (any(is.na(info[,])))
    stop("Failed to compute x-coordinates.")
  if (!silent) print("Computing y-coordinates...")
  y_4col= cell_yavg(1:cat$nrows,cat$nrows,cat$yllcorner,cat$cellsize)
  ycoords= matrix(y_4col, ncol=cat$ncols, nrow=cat$nrows, byrow=FALSE)
  tmp= tapply(X=ycoords, INDEX=cat$matrix, FUN=mean)
  tmp= data.frame(id=names(tmp), y=round(as.numeric(tmp)))
  tmp= subset(tmp, tmp$id != cat$nodata_value)
  info= merge(x= info, y=tmp, by.x="object", by.y="id", all=TRUE)
  rm(y_4col, ycoords, tmp)
  if (any(is.na(info[,])))
    stop("Failed to compute y-coordinates.")
  # Compute area
  if (!silent) print("Computing areas...")
  tmp= tapply(X=cat$matrix, INDEX=cat$matrix, FUN=length) * (cat$cellsize^2)
  tmp= data.frame(id=names(tmp), area=round(as.numeric(tmp)))
  tmp= subset(tmp, tmp$id != cat$nodata_value)
  info= merge(x= info, y=tmp, by.x="object", by.y="id", all=TRUE)
  rm(tmp)
  if (any(is.na(info[,])))
    stop("Failed to compute catchment areas.")
  # Add prefix to IDs
  info$object= paste(prefix, info$object, sep="")
  # Return data frame
  return(info)
}

