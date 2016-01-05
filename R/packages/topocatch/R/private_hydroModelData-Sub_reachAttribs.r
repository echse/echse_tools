
# This function computes reach attributes.
# It is called by the pulic high-level function 'hydroData'.

reachAttribs= function(
  tab_objDecl,      # object declaration table (output of 'hydrolinkage')
  tab_downObj,      # table listing direct downstream neighbors (output of 'hydrolinkage')
  tab_catchAttribs, # table with catchment attributes (output of 'catchAttribs')
  fileDEM,          # Grid file with DEM
  fileSHP,          # Shape file with drainage lines
  id_field,         # ID field in the shape file
  class_field,      # Class field in the shape file
  classname_reach,     # Class name for reach objects
  classname_minireach, # Class name for minireach objects
  updateSHP,           # Add fields to shape file's attr. table?
  min_slope= 0.0001,   # slope to be used where computed value is zero
  silent=TRUE          # Switch for diagnostic messages
) {
  # Check args
  checkArg(arg=tab_objDecl, len=NULL, type="data.frame")
  if (!all(c("object","objectGroup") %in% names(tab_objDecl)))
    stop("Missing columns in object declaration table.")
  checkArg(arg=tab_downObj, len=NULL, type="data.frame")
  if (!all(c("id","id_down") %in% names(tab_downObj)))
    stop("Missing columns in table of downstream neighbors.")
  checkArg(arg=tab_catchAttribs, len=NULL, type="data.frame")
  if (!all(c("object","area") %in% names(tab_catchAttribs)))
    stop("Missing columns in catchment attributes table.")
  checkFileIn(fileDEM)
  checkFileIn(fileSHP)
  checkArg(arg=id_field, len=1, type="character")
  checkArg(arg=class_field, len=1, type="character")
  checkArg(arg=classname_reach, len=1, type="character")
  checkArg(arg=classname_minireach, len=1, type="character")
  checkArg(arg=updateSHP, len=1, type="logical")
  checkArg(arg=min_slope, len=1, type="numeric")
  checkArg(arg=silent, len=1, type="logical")

  # Init reach attributes table (don't use the shape file but the declaration table
  # because head reaches may have been deleted during linkage)
  inds= which(tab_objDecl$objectGroup %in% c(classname_reach,classname_minireach))
  info= data.frame(object= tab_objDecl$object[inds], class= tab_objDecl$objectGroup[inds])

  # Read full line data
  if (!silent) print("Converting shape file...")
  tab_shp= lineShapeAsTable(file=fileSHP, id_field=id_field, attribs=c(class_field), endsOnly=FALSE)
  tab_shp= merge(x=tab_shp$shp, y=tab_shp$att, by=id_field, all=TRUE)
  if (any(is.na(tab_shp)))
    stop("Failed to merge feature geometry and attributes.")

  if (!silent) print("Computing length of reach objects...")
  tab_shp$xy= paste(tab_shp$x, tab_shp$y, sep="_")
  tmp= tapply(X=tab_shp$xy, INDEX=tab_shp$id, FUN=lineLength, sepchar="_")
  tmp= data.frame(id=names(tmp), length=round(as.numeric(tmp),3))
  info= merge(x=info, y=tmp, by.x="object", by.y="id", all.x=TRUE, all.y=FALSE)
  rm(tmp)
  if (any(is.na(info[,]))) stop("Failed compute length of reach objects.")

  # Re-read reach data (now only endpoints of lines)
  rm(tab_shp)
  tab_shp= lineShapeAsTable(file=fileSHP, id_field=id_field, attribs=c(class_field), endsOnly=TRUE)
  tab_shp= merge(x=tab_shp$shp, y=tab_shp$att, by=id_field, all=TRUE)
  if (any(is.na(tab_shp)))
    stop("Failed to merge feature geometry and attributes.")

  if (!silent) print("Saving coordinates of reach objects...")
  info= merge(x=info, y=tab_shp[,c("id","x1","y1","x2","y2")], by.x="object", by.y="id", all.x=TRUE, all.y=FALSE)
  if (any(is.na(info[,])))
    stop("Failed to save coordinates of reach objects.")
  for (i in c("x1","y1","x2","y2")) {
    info[,i]= round(info[,i],3)
  }

  if (!silent) print("Computing elevations of reach objects...")
  dem= geogrid.readAscii(fileDEM)
  # Determine elevation of end points
  z1= geogrid.valuesAtPoints(grid=dem, xy_table=data.frame(x=tab_shp[,"x1"], y=tab_shp[,"y1"]))
  z2= geogrid.valuesAtPoints(grid=dem, xy_table=data.frame(x=tab_shp[,"x2"], y=tab_shp[,"y2"]))
  if (any(z1 == dem$nodata_value) || any(z2 == dem$nodata_value))
    stop("Failed to set elevations for reach end(s). Reach(es) not within valid DEM area.")
  tab_shp= cbind(tab_shp, elev_min= pmin(z1, z2))
  tab_shp= cbind(tab_shp, elev_max= pmax(z1, z2))
  # Add elevation range to attributes table
  info= merge(x=info, y=tab_shp[,c("id","elev_min","elev_max")], by.x="object", by.y="id", all.x=TRUE, all.y=FALSE)
  if (any(is.na(info[,])))
    stop("Failed to compute elevations of reach objects.")
  # Clean up
  rm(dem)
  rm(z1)
  rm(z2)
  rm(tab_shp)
  # Add slope column
  if (!silent) print("Computing slope of reach objects...")
  info= cbind(info, slope= (info$elev_max - info$elev_min) / info$length)
  info$slope[info$slope <= 0]= min_slope
  info$slope= signif(info$slope, 3)

################################################################################

  if (!silent) print("Computing upstream catchment of reach objects...")
  # Identify name of object downstream (I.E. OUTSIDE) of the system
  ix= which(!(tab_downObj$id_down %in% tab_downObj$id))
  if (length(ix) != 1) 
    stop("Failed to identify object downstream of system.")

  ix= ix[1]
  # Add artificial end object having itself as downstream neighbor (circle)
  ID_END= tab_downObj$id_down[ix]
  tab_downObj= rbind(tab_downObj, list(id=ID_END, id_down=ID_END))

  # Initialize upstream area for catchment objects
  tab_downObj= merge(x=tab_downObj, y=tab_catchAttribs[,c("object","area")], by.x="id", by.y="object", all=TRUE)
  # Initialize upstream area to zero for all other objects
  tab_downObj$area[is.na(tab_downObj$area)]= 0
  if (any(is.na(tab_downObj)))
    stop("Failed to add area field to table of downstream neighbors.")

  # Determine index of artificial end object
  IX_END= match(ID_END, tab_downObj$id)

  # Compute areas
  tab_downObj$area_send= tab_downObj$area
  while (sum(tab_downObj$area_send) > 0) {
    tab_downObj$area_reci= 0
    tmp= tapply(X=tab_downObj$area_send, INDEX=tab_downObj$id_down, FUN=sum)
    inds_reci= match(names(tmp), tab_downObj$id)
    tab_downObj$area_reci[inds_reci]= as.numeric(tmp)
    tab_downObj$area= tab_downObj$area + tab_downObj$area_reci
    tab_downObj$area_send= tab_downObj$area_reci
    tab_downObj$area_send[IX_END]= 0
  }
  # Remove/rename columns
  tab_downObj$id_down= NULL
  tab_downObj$area_send= NULL
  tab_downObj$area_reci= NULL
  names(tab_downObj)[names(tab_downObj) == "area"]= "upstr_area"
  # Adjust units  
  tab_downObj$upstr_area= tab_downObj$upstr_area / 1.e06

  # Merge with attributes table
  info= merge(x=info, y=tab_downObj, by.x="object", by.y="id", all.x=TRUE, all.y=FALSE)
  if (any(is.na(info[,])))
    stop("Failed to compute upstream catchment of reach objects.")

################################################################################

  # Update the shape file's attribute table (we can add info for non-head reaches only)
  if (updateSHP) {
    if (!silent) print("Updating shape file's attribute table...")
    # Read attr. table
    dbfname= paste(substr(fileSHP,1,nchar(fileSHP)-3),"dbf",sep="")
    attTab= shapefiles::read.dbf(dbfname, header=TRUE)
    attTab= attTab$dbf
    # Delete fields that will be updated if they already exist (may happen in repeated calls, for example)
    del= which((names(attTab) %in% names(info)) & (!(names(attTab) %in% c(id_field,class_field))))
    attTab[,del]=NULL
    # Add reach attributes (but avoid duplicated class field)
    newTab= merge(x=attTab, y=info[,names(info) != "class"], by.x=id_field, by.y="object", all=TRUE)
    # Set nodata value for non-reach objects and those that are not modelled (dropped from system due to missing upstream input)
    newTab[is.na(newTab)]= -9999
    # Check and write result
    if (nrow(newTab) != nrow(attTab))
      stop("Failed to update shape file's attribute table.")
    foreign::write.dbf(newTab, dbfname)   # The 'foreign' package is loaded via 'shapefiles' (but this function is masked, therefore we use ::)
  }
  return(info)
}

