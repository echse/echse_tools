
#' Pre-processing of geo data for hydrological catchment modeling
#'
#' Type \code{help(package="topocatch")} to inspect the package description.
#'
#' @name topocatch-package
#' @aliases topocatch
#' @docType package
{}

################################################################################

#' Filling of sinks in a digital elevation model (DEM)
#'
#' See the workhorse function \code{\link{sinkfill}} for details.
#'
#' @param fileIn Name/path of input file (ASCII grid).
#' @param fileOut Name/path of output file (ASCII grid).
#' @inheritParams sinkfill
#' @inheritParams geogrid.readAscii
#' @inheritParams geogrid.writeAscii
#'
#' @return \code{NULL}
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

dem.fill= function(fileIn, fileOut, ndigits, replace=FALSE, silent=TRUE
) {
  # Check args
  # Bool
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")
  # Files
  checkFileIn(fileIn)
  checkFileOut(fileOut, replace=replace)
  # Numbers
  checkArg(arg=ndigits, len=1, type="integer")
  # Process
  if (!silent) print("Reading DEM grid...")
  dem= geogrid.readAscii(fileIn)
  if (!silent) print("Filling sinks...")
  dem2= sinkfill(grid=dem, ndigits=ndigits, silent=silent)
  if (!silent) print("Writing sink-filled DEM grid...")
  geogrid.writeAscii(grid=dem2, file=fileOut, replace=replace)
  return(invisible(NULL))
}

################################################################################

#' Analyze a digital elevation model (DEM)
#'
#' See the workhorse functions \code{\link{flowdir}} and
#' \code{\link{concTimeIndex}} for details.
#'
#' @param fileDEM Name/path of the INPUT file containing the \emph{sink-filled}
#'   DEM (ASCII grid).
#' @param fileDIR Name/path of the OUTPUT file containing flow direction codes
#'   (ASCII grid).
#' @param fileACC Name/path of the OUTPUT file containing flow accumulation data
#'   (ASCII grid).
#' @param fileCTI Name/path of the OUTPUT file containing values of the
#'   concentation time index (ASCII grid).
#' @param fileSHP Name/path of the OUTPUT file containing the generated river
#'   net (shape file format).
#' @param minlength_reach Minimum length of a river section to be classified as
#'   a reach. Shorter sections are considered as so-called 'mini-reaches'.
#' @param classname_reach Class name for reach objects.
#' @param classname_minireach Class name for mini-reach objects (see the
#'   \code{minlength_reach} argument).
#' @param id_field Name of the ID field in the output shape file.
#' @param class_field Name of the class field in the output shape file.
#' @inheritParams flowdir
#' @inheritParams concTimeIndex
#' @inheritParams flowPaths
#' @inheritParams geogrid.readAscii
#' @inheritParams geogrid.writeAscii
#'
#' @return ID of the system's outlet, i.e. the most downstream reach (integer).
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

dem.analyze= function(fileDEM, fileDIR, fileACC, fileCTI, fileSHP,
  crit_source_area,
  x_inBasin=NULL, y_inBasin=NULL, id_field="id", class_field="class",
  minlength_reach=100, classname_reach="rch", classname_minireach="minirch",
  dz_min=0.1, replace=FALSE, silent=TRUE
) {
  # Check args
  # Bool
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")
  # Files
  checkFileIn(fileDEM)
  checkFileOut(fileDIR, replace=replace)
  checkFileOut(fileACC, replace=replace)
  checkFileOut(fileCTI, replace=replace)
  checkFileOut(fileSHP, replace=replace)
  # Names
  checkArg(arg=id_field, len=1, type="character")
  checkArg(arg=class_field, len=1, type="character")
  checkArg(arg=classname_reach, len=1, type="character")
  checkArg(arg=classname_minireach, len=1, type="character")
  # Numbers
  checkArg(arg=crit_source_area, len=1, type="numeric")
  checkArg(arg=minlength_reach, len=1, type="numeric")
  checkArg(arg=dz_min, len=1, type="numeric")
  # Process
  if (!silent) print("Reading DEM grid...")
  dem= geogrid.readAscii(fileDEM)
  if (!silent) print("Computing flow direction codes...")
  fdir= flowdir(grid=dem, silent=silent)
  if (!silent) print("Writing grid of flow direction codes...")
  geogrid.writeAscii(grid=fdir, file=fileDIR, replace=replace)
  if (!silent) print("Computing flow accumulation...")
  facc= flowacc(grid=fdir, silent=silent)
  if (!silent) print("Writing grid of flow accumulation...")
  geogrid.writeAscii(grid=facc, file=fileACC, replace=replace)
  if (!silent) print("Computing concentration time indices...")
  cti= concTimeIndex(grid_dem=dem, grid_flowdir=fdir, grid_flowacc=facc,
    crit_source_area=crit_source_area, dz_min=dz_min, silent=silent)
  if (!silent) print("Writing grid of concentration time indices...")
  cti$matrix= round(cti$matrix,1)
  geogrid.writeAscii(grid=cti, file=fileCTI, replace=replace)
  rm(cti) # Clean up
  if (!silent) print("Computing flow paths...")
  tmp= flowPaths(grid_flowdir=fdir, grid_flowacc=facc, crit_source_area=crit_source_area,
    x_inBasin=x_inBasin, y_inBasin=y_inBasin, silent=silent)
  id_outlet= tmp$id_outlet # Set outlet ID
  rm(fdir) # Clean up
  rm(facc) # Clean up
  if (!silent) print("Assembling attribute table...")
  tmp= tmp$shpTable # Only keep the table in tmp
  tmp$xy= paste(tmp$x, tmp$y, sep="_")
  len= tapply(X=tmp$xy, INDEX=tmp$id, FUN=lineLength, sepchar="_")
  tmp$xy= NULL
  attr= data.frame(id=as.integer(names(len)), len=as.numeric(len), class=NA)
  attr$class[attr$len >= minlength_reach]= classname_reach
  attr$class[attr$len < minlength_reach]= classname_minireach
  attr$len= NULL
  names(attr)[which(names(attr) == "id")]= id_field
  names(attr)[which(names(attr) == "class")]= class_field
  if (!silent) print(paste("Identified ",sum(attr[,class_field]==classname_reach),
    " object(s) of class '",classname_reach,"'",sep=""))
  if (!silent) print(paste("Identified ",sum(attr[,class_field]==classname_minireach),
    " object(s) of class '",classname_minireach,"'",sep=""))
  if (!silent) print("Assembling shape file...")
  names(tmp)[which(names(tmp)== "id")]= id_field
  s= convert.to.shapefile(shpTable=tmp, attTable=attr, field=id_field, type=3)
  if (nchar(fileSHP) > 4) {
    if (substr(fileSHP,nchar(fileSHP)-3,nchar(fileSHP)) == ".shp")
      fileSHP= substr(fileSHP,1,nchar(fileSHP)-4)
  }
  if (file.exists(paste(fileSHP,".shp",sep="")) && (!replace))
    stop(paste("Shape file '",paste(fileSHP,".shp",sep=""),"' already exists.",sep=""))
  if (!silent) print("Writing shape file...")
  write.shapefile(shapefile=s, out.name=fileSHP, arcgis=TRUE)
  # Return outlet ID
  return(id_outlet)
}

################################################################################

#' Derive input for hydrological catchment modeling from pre-processed geo data
#'
#' The function identifies the relevant objects for object-based hydrological
#' catchment modeling from pre-processed geo data sets. It creates several output
#' files which can directly serve as an input for hydrological models, namely those
#' build with the ECHSE simulation environment.
#'
#' @param fileSHP Name/path of an INPUT shape file with line features representing
#'   drainage lines. This file is either an output of \code{\link{dem.analyze}},
#'   a manually edited version of such output or an external file (possibly
#'   created by manual digitizing). If \code{updateSHP} is \code{TRUE}, some new
#'   fields will be appended to the shape file's attribute table (making it an
#'   input \emph{and} output).
#' @param fileDEM Name/path of the INPUT file containing the \emph{sink-filled}
#'   DEM (ASCII grid). It is used here to estimate the bed slope of river reaches.
#'   If the shape file supplied as \code{fileSHP} was created by a call to
#'   \code{\link{dem.analyze}}, the same DEM as in this call should be used.
#' @param fileDIR Name/path of the INPUT file containing flow direction codes
#'   computed by \code{\link{dem.analyze}} (ASCII grid).
#' @param fileCAT Name/path of the OUTPUT grid file showing shape, position and
#'   extent of the created catchments (ASCII grid). Each catchment is identified
#'   by an ID number. This integer code is derived from the IDs of the
#'   corresponding features in the input shape file (\code{fileSHP}).
#' @param fileAttrCAT Name/path of a tabular OUTPUT file listing basic attributes
#'   of the generated catchments such as the areal extent (in field 'area') and the
#'   positions of the center of gravity ('x' and 'y' fields). The unit of the area
#'   is square of the basic length unit used in the input files. See also the
#'   \code{findAffectedGages} argument.
#' @param fileAttrRCH Name/path of a tabular OUTPUT file listing basic attributes
#'   of reach objects (incl. 'mini-reaches'). The attributes include the reach
#'   length (field 'length'), the elevation of the two reach ends (fields
#'   'elev_min' and 'elev_max'), a rough estimate of the bed slope derived from
#'   the former information ('slope' field) as well as the total area of the
#'   reach's upstream catchment (in field 'upstreamArea'). See also the
#'   \code{findAffectedGages} argument.
#' @param fileObjDecl Name/path of a tabular OUTPUT file to be used as an input
#'   by ECHSE-based hydrological models. This file represents the so-called
#'   'object declaration table' holding IDs and class info for all objects.
#' @param fileObjLink Name/path of a tabular OUTPUT file to be used as an input
#'   by ECHSE-based hydrological models. This file represents the so-called
#'   'object linkage table'. It provides info on object interactions.
#' @param id_outlet ID of the line feature representing the system's outlet
#'   (ID of the most downstream reach). The value must exist in the ID field of
#'   the shape file's attribute table specified by \code{id_field}. If the shape
#'   file supplied as \code{fileSHP} was created by a call to \code{\link{dem.analyze}},
#'   without subsequent modification, the return value of this function is an
#'   appropriate input.
#' @param id_field Name of the field in the shape file's attribute table
#'   containing feature IDs.
#' @param class_field Name of the field in the shape file's attribute table
#'   holding information on the the features' classes. This field is used, for
#'   example, to distinguish between reaches, pipes, reservoirs, etc. See the
#'   \code{classes_with_catchment} argument for further details.
#' @param classname_reach Class name used for reach objects. If the shape
#'   file supplied as \code{fileSHP} was created by a call to
#'   \code{\link{dem.analyze}}, the value should be consistent with the one used
#'   in that call.
#' @param classname_minireach Class name used for very short reach objects with
#'   negligible travel time.  If the shape file supplied as \code{fileSHP} was
#'   created by a call to \code{\link{dem.analyze}}, the value should be
#'   consistent with the one used in that call.
#' @param classname_node Class name to be assigned to node objects (junctions).
#' @param classname_catch Class name to be assigned to catchment objects.
#' @param classname_gage Class name to be assigned to gage objects, if existent.
#' @param classes_with_catchment A vector of class names (strings). Catchments are
#'   generated only for features belonging to those classes. These class names are
#'   expected to exist in the attribute table's field specified by \code{class_field}.
#'   Typically, catchments are to be generated for objects of the reach class, at
#'   least. Lake and reservoir objects, usually have a catchment on its own as well.
#' @param nbuffer An integer value to control the thickness of the rasterized
#'   lines during vector-to-raster conversion. The converted lines will be
#'   \eqn{nbuffer * 2 + 1} cells wide. Thus, the default value of 1 results in
#'   rasterized lines being 3 cells wide.
#' @param coord_tol A small distance value (same spatial unit as used in the
#'   shape file) to account for precision problems. The positions of two lines'
#'   end points are assumed to be identical if their distance is <= the value of
#'   \code{coord_tol}. The default is 1.e-03.
#'   Shape files created by a call to \code{\link{dem.analyze}} are always
#'   \emph{exact} in this respect and the small default value is appropriate.
#'   If the shape file was created manually in a GIS without the use of a snapping
#'   mechanism (which is \emph{not} recommended), a larger value might be required.
#'   The value of \code{coord_tol} is also used when gage objects are specified by
#'   position via the \code{gageLocations} argument.
#' @param prefix_node String used as a prefix when generating the IDs of node
#'   objects based on the features ID of the respective downstream objects.
#' @param prefix_catch String used as a prefix when generating the IDs of
#'   catchment objects based on the IDs of the corresponding features from the
#'   shape file.
#' @param min_slope Minimum slope of river bed (dimensionless). This value is
#'   used as a substitute for zero-slope values in the reach attributes table.
#'   Such zero-values typically exist in nearly flat areas.
#' @param updateSHP If \code{TRUE}, fields will be appended to the shape file's
#'   attribute table. Those fields contain information for reach-like objects
#'   such as the length and the total upstream area. The value of -9999 is used
#'   in these fields for non-reach-like objects and those objects which were
#'   dropped from the system because of missing inflow from upstream.
#' @param namesIO A data frame with 2 columns named 'target' and 'source'. Each
#'   record represents the names of a pair of corresponding input-output variables. A typical
#'   example would be a data frame with just a single record and the strings
#'   "inflow" and "outflow" in the 'target' and 'source' column, respectively.
#'   Additional records in the data frame declare additional I/O variables whose
#'   values are exchanged between interacting objects.
#' @param gageLocations A data frame with (at least) 3 columns named 'id',
#'   'x' and 'y'. Each record defines the position of a gage. By default, this
#'   data frame is empty. The positions given by x and y \emph{must} coindide
#'   with the coordinates of the start or end point of a line feature in the
#'   shape file. Thus, a gage can only be located at either end of a reach but not
#'   somewhere in the mid. Note that the coordinates must match to the precision
#'   defined by \code{coord_tol}.
#'  @param findAffectedGages If \code{TRUE}, the function identifies the gage
#'    objects being affected by each particular object. For each gage, an additional
#'    column is appended to the output files \code{fileObjDecl}, \code{fileAttrCAT},
#'    and \code{fileAttrRCH}. Column names are taken from the 'id'
#'    field of \code{gageLocations}. The values in those columns are 1 or 0, where
#'    1 means that the gage (column) is affected by the simulated object (row) and
#'    0 indicates the opposite (no interaction). The default for \code{findAffectedGages}
#'    is \code{FALSE} because the current algorithm may be slow for very large systems.
#'    If \code{gageLocations} is empty this argument setting has no effect.
#'  @param replace Should existing output files be silently replaced?
#'  @param silent Logical value to turn diagnostic messages on/off.
#'
#' @return The function returns \code{NULL}. All computed information is written
#'   to the respective output files.
#'
#' @note Not all of the features in the shape file are retained as reach-like
#'   objects. Those features that do not have a catchment upstream are dropped
#'   automatically during creation of the various output tables. This is true in
#'   particular for head-reaches.
#'
#'   If gage objects are to be considered, an iterative approach may be convienent.
#'   In the 1st step, this function is called without any gage specifications (default for
#'   \code{gageLocations}). Suitable gage positions can then be identified in the GIS
#'   based on the updated shape file if \code{updateSHP} was set to \code{TRUE}
#'   (coordinates are available as attributes and can be queried). In the 2nd step,
#'   the function is called with gage positions supplied in \code{gageLocations}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export

hydroModelData= function(
  fileSHP,
  fileDEM,
  fileDIR,
  fileCAT,
  fileAttrCAT,
  fileAttrRCH,
  fileObjDecl,
  fileObjLink,
  id_outlet,
  id_field="id",
  class_field="class",
  classname_reach="rch",
  classname_minireach="minirch",
  classname_node="node",
  classname_catch="cat",
  classname_gage="gage",
  classes_with_catchment= c(classname_reach),
  nbuffer=1,
  coord_tol=1.e-03,
  prefix_node="node_",
  prefix_catch="cat_",
  min_slope= 0.0001,
  updateSHP= FALSE,
  namesIO= data.frame(target=c("qi"),source=c("qx")),
  gageLocations= data.frame(id=c(),x=c(),y=c()),
  findAffectedGages= FALSE,
  replace=FALSE,
  silent=TRUE
){
  # Check args
  checkFileIn(fileSHP)
  checkFileIn(fileDEM)
  checkFileIn(fileDIR)
  checkFileOut(fileCAT, replace=replace)
  checkFileOut(fileAttrCAT, replace=replace)
  checkFileOut(fileAttrRCH, replace=replace)
  checkFileOut(fileObjDecl, replace=replace)
  checkFileOut(fileObjLink, replace=replace)
  checkArg(arg=id_outlet, len=1, type="integer")
  checkArg(arg=id_field, len=1, type="character")
  checkArg(arg=class_field, len=1, type="character")
  checkArg(arg=classname_reach, len=1, type="character")
  checkArg(arg=classname_minireach, len=1, type="character")
  checkArg(arg=classname_node, len=1, type="character")
  checkArg(arg=classname_catch, len=1, type="character")
  checkArg(arg=classname_gage, len=1, type="character")
  checkArg(arg=classes_with_catchment, len=NULL, type="character")
  checkArg(arg=nbuffer, len=1, type="integer")
  checkArg(arg=coord_tol, len=1, type="numeric")
  checkArg(arg=prefix_node, len=1, type="character")
  checkArg(arg=prefix_catch, len=1, type="character")
  checkArg(arg=min_slope, len=1, type="numeric")
  checkArg(arg=updateSHP, len=1, type="logical")
  checkArg(arg=namesIO, len=NULL, type="data.frame")
  checkArg(arg=gageLocations, len=NULL, type="data.frame")
  checkArg(arg=findAffectedGages, len=1, type="logical")
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")

  # Create catchment grid
  if (!silent) print("Computing catchments...")
  catchments(fileDIR=fileDIR, fileSHP=fileSHP, fileCAT=fileCAT,
    id_field=id_field, class_field=class_field,
    classes_with_catchment= classes_with_catchment,
    nbuffer=nbuffer, replace=replace, silent=silent
  )
  # Create data frame with basic catchment attributes  
  if (!silent) print("Computing catchment attributes...")
  tab_catchAttribs= catchAttribs(
    fileCAT=fileCAT, prefix=prefix_catch, silent=silent
  )
  # Analyze drainage network
  if (!silent) print("Analyzing drainage network...")
  info= hydroLinkage(shapefile=fileSHP, id_field=id_field, class_field=class_field,
    coord_tol=coord_tol, id_outlet=id_outlet,
    classes_with_catchment=classes_with_catchment,
    classname_reach=classname_reach, classname_minireach=classname_minireach,
    classname_node=classname_node, classname_catch=classname_catch, classname_gage=classname_gage,
    prefix_node=prefix_node, prefix_catch=prefix_catch,
    tab_var=namesIO, tab_gages=gageLocations, silent=silent
  )
  # Create data frame with reach attributes
  if (!silent) print("Computing reach attributes...")
  tab_reachAttribs= reachAttribs(
    tab_objDecl=info$objDecl,
    tab_downObj=info$downObj,
    tab_catchAttribs=tab_catchAttribs,
    fileDEM=fileDEM, fileSHP=fileSHP, id_field=id_field, class_field=class_field,
    classname_reach=classname_reach, classname_minireach=classname_minireach,
    updateSHP=updateSHP, min_slope=min_slope, silent=silent
  )
  if (findAffectedGages) {
    colname_affectedGagesCode= "code_affectedGages"
    if (!silent) print("Identifying affected gages...")
    # Identify the gages affected by each object
    if (nrow(gageLocations) > 0) {
      # Force string conversion
      gageLocations$id= as.character(gageLocations$id)
      # Check for name conflicts
      reserved= unique(c(id_field,class_field,
        names(info$objDecl), names(tab_catchAttribs), names(tab_reachAttribs)))
      if (any(gageLocations$id %in% reserved))
        stop(paste("Name conflict. A gage name must not be one of '",
          paste(reserved,collapse="', '"),"'.",sep=""))
      if (colname_affectedGagesCode %in% reserved)
        stop(paste("Error in package's source code. Column '",colname_affectedGagesCode,
          "' already exists.",sep=""))
      # Update the object declaration table
      info$objDecl= assignGages(tab_objDecl= info$objDecl, tab_objLink= info$objLink,
        gageObjects= gageLocations$id)
      # Add the info also merged into a single field. Let it be a character field to avoid
      # problems with the representation of very large numbers (in the case of many gages)
      info$objDecl[,colname_affectedGagesCode]= apply(X=info$objDecl[,gageLocations$id],
        MARGIN=1,FUN=paste,collapse="")
      info$objDecl[,colname_affectedGagesCode]= paste("code_",info$objDecl[,colname_affectedGagesCode],sep="")
      # Update the catchment attribute table
      tab_catchAttribs= merge(x=tab_catchAttribs, y=info$objDecl[,c("object",gageLocations$id,colname_affectedGagesCode)],
        by="object", all.x=TRUE, all.y=FALSE)
      if (any(is.na(tab_catchAttribs)))
        stop("Failed to add gage info to catchment attribute table.")
      # Update the reach attribute table
      tab_reachAttribs= merge(x=tab_reachAttribs, y=info$objDecl[,c("object",gageLocations$id,colname_affectedGagesCode)],
        by="object", all.x=TRUE, all.y=FALSE)
      if (any(is.na(tab_reachAttribs)))
        stop("Failed to add gage info to reach attribute table.")
      # Update the shape file's attribute table if requested
      if (updateSHP) {
        if (!silent) print("Updating shape file's attribute table...")
        # Read attr. table
        dbfname= paste(substr(fileSHP,1,nchar(fileSHP)-3),"dbf",sep="")
        attTab= read.dbf(dbfname, header=TRUE)
        attTab= attTab$dbf
        # Delete fields that will be updated if they already exist (may happen in repeated calls, for example)
        del= which(names(attTab) %in% c(gageLocations$id,colname_affectedGagesCode))
        attTab[,del]=NULL
        # Add gage info
        newTab= merge(x=attTab, y=info$objDecl[,c("object",gageLocations$id,colname_affectedGagesCode)], by.x=id_field,
          by.y="object", all.x=TRUE, all.y=FALSE)
        # Set info for non-simulated objects to special value
        for (i in 1:nrow(gageLocations)) {
          newTab[is.na(newTab[,gageLocations$id[i]]),gageLocations$id[i]]= -9999
        }
        newTab[is.na(newTab[,colname_affectedGagesCode]),colname_affectedGagesCode]= -9999
        foreign::write.dbf(newTab, dbfname)   # The 'foreign' package is loaded via 'shapefiles' (but this function is masked, therefore we use ::)
      }
    }
  }
  # Write output tables  
  write.table(x=tab_catchAttribs, file=fileAttrCAT, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  write.table(x=tab_reachAttribs, file=fileAttrRCH, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  write.table(x=info$objDecl, file=fileObjDecl, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  write.table(x=info$objLink, file=fileObjLink, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  # Return nothing
  return(invisible(NULL))
}


