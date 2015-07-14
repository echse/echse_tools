
# This function identifies the structure of a river network.
# It is called by the pulic high-level function 'hydroData'.

hydroLinkage = function (
  shapefile,                 # shape file with drainage lines
  id_field,                  # id field in shape file
  class_field,               # class field in shape file
  coord_tol,                 # tolerance when checking for identical spatial position
  id_outlet,                 # ID of the feature representing the catchment's outlet
  classes_with_catchment,    # classes to which a catchment should be assigned
  classname_reach,           # class name for reach objects
  classname_minireach,       # class name for minireach objects
  classname_node,            # basic class name for node objects
  classname_catch,           # basic class name for catchment objects
  classname_gage,            # class name for gage objects
  prefix_node,               # prefix to be used in IDs of node objects
  prefix_catch,              # prefix to be used in IDs of catchment objects
  tab_var,                   # table of names for the exchanged variables
  tab_gages,                 # table with gage-coordinates
  silent=TRUE                # switch for diagnostic messages
) {

  # Define helper function
  uplist= function(id_this, vect_idDown, vect_id) {
    vect_id[vect_idDown %in% id_this]
  }

  # Define constants
  SEPCHAR=";"               # character to separate IDs in a string-list
  COORDSEP="_"              # separator in position string between x y

  #.............................................................................
  # Check arguments
  #.............................................................................

  # Basic checks
  checkFileIn(shapefile)
  checkArg(arg=id_field,len=1,type="character")
  checkArg(arg=class_field,len=1,type="character")
  checkArg(arg=coord_tol,len=1,type="numeric")
  checkArg(arg=id_outlet,len=1,type="integer")
  checkArg(arg=classes_with_catchment,len=NULL,type="character")
  checkArg(arg=classname_reach,len=1,type="character")
  checkArg(arg=classname_minireach,len=1,type="character")
  checkArg(arg=classname_node,len=1,type="character")
  checkArg(arg=classname_catch,len=1,type="character")
  checkArg(arg=classname_gage,len=1,type="character")
  checkArg(arg=prefix_node,len=1,type="character")
  checkArg(arg=prefix_catch,len=1,type="character")
  checkArg(arg=tab_var,len=NULL,type="data.frame")
  checkArg(arg=tab_gages,len=NULL,type="data.frame")
  checkArg(arg=silent,len=1,type="logical")

  # More checks: Name prefixes
  if (prefix_catch == prefix_node)
     stop("Arguments 'prefix_catch' and 'prefix_node' must not have the same value.")

  # More checks: Table of variable names
  cols= c("source","target")
  if (!all(cols %in% names(tab_var)))
    stop(paste("Argument 'tab_var' must be a data frame with columns '",
      paste(cols,collapse="', '"),"'.",sep=""))

  # More checks: Table of gages
  if (nrow(tab_gages) > 0) {
    cols= c("id","x","y")
    if (!all(cols %in% names(tab_gages)))
      stop(paste("Argument 'tab_gages' must be a data frame with columns '",
        paste(cols,collapse="', '"),"'.",sep=""))
    if (any(duplicated(tab_gages[,c("x","y")])))
      stop("Duplicate gage positions in 'tab_gages'.")
    if (any(duplicated(tab_gages$id)))
      stop("Duplicate gage IDs in 'tab_gages'.")
  }

  #.............................................................................
  # Convert shapefile to table (line endings only)
  #.............................................................................

  if (!silent) print("Converting shapefile...")
  tab_net= lineShapeAsTable(file=shapefile, id_field=id_field, attribs=c(class_field), endsOnly=TRUE)
  tab_net= merge(x=tab_net$shp, y=tab_net$att, by=id_field, all=TRUE)
  if (any(is.na(tab_net)))
    stop("Failed to merge feature geometry and attributes.")
  # Rename class field
  names(tab_net)[which(names(tab_net) == class_field)]= "class"

  #.............................................................................
  # Initial checks and conversions
  #.............................................................................

  if (!silent) print("Checking integrity of input data...")

  # Convert ID column to character
  tab_net$id= as.character(tab_net$id)
  # Check for unique IDs
  if (length(unique(tab_net$id)) != nrow(tab_net))
    stop("Duplicates detected in ID column. Values must be unique.")

  # Round coordinates according to the specified tolerance
  tab_net$x1= round(tab_net$x1 / coord_tol)
  tab_net$y1= round(tab_net$y1 / coord_tol)
  tab_net$x2= round(tab_net$x2 / coord_tol)
  tab_net$y2= round(tab_net$y2 / coord_tol)
 
  # Merge x and y info in single columns (facilitates comparison of locations)
  tab_net$pos1= paste(as.character(tab_net$x1), as.character(tab_net$y1), sep=COORDSEP)
  tab_net$pos2= paste(as.character(tab_net$x2), as.character(tab_net$y2), sep=COORDSEP)
  tab_net$x1= NULL
  tab_net$y1= NULL
  tab_net$x2= NULL
  tab_net$y2= NULL

  # Set character to separate element names in a list. This character must not
  # be part of an element's ID
  if (length(grep(SEPCHAR,tab_net$id,fixed=TRUE)) != 0) {
    stop("Character '",SEPCHAR,"' must not be present in element IDs.")
  }

  # Make sure that class names for catchments and nodes are not in use already
  if (any(tab_net$class == classname_node))
    stop(paste("Input already contains elements of class '",classname_node,"'.",sep=""))
  if (any(tab_net$class == classname_catch))
    stop(paste("Input already contains elements of class '",classname_catch,"'.",sep=""))

  # Make sure that the classes with corresponding catchments exist
  for (i in 1:length(classes_with_catchment)) {
    if (!(classes_with_catchment[i] %in% unique(tab_net$class)))
    stop(paste("Assignment of catchments to elements of class '",classes_with_catchment[i],
      "' is requested but this class is not present in the input table.",sep=""))
  }

  # Make sure that a catchment is assigned to objects of the reach class
  # and that NO catchment is assigned to objects of the minireach class
  if (!(classname_reach %in% classes_with_catchment))
    stop("Name of the reach class is missing in the list of classes with a catchment.")
  if (classname_minireach %in% classes_with_catchment)
    stop("Name of the minireach class must not appear in the list of classes with a catchment.")

  #.............................................................................
  # Identify the outlet (last line in the system) and perform checks/corrections
  #.............................................................................

  if (!silent) print("Identifying the system's outlet...")

  # Find element by ID
  ix_outlet= which(tab_net$id == id_outlet)
  if (length(ix_outlet) != 1) {
    stop(paste("Found ",length(ix_outlet)," potential system outlet(s) with",
      " ID '",id_outlet,"'.",sep=""))
  }
  # One end of the outlet must be linked and the other end must NOT be linked
  inds= which(tab_net$id != id_outlet)
  linked_1= any( (tab_net$pos1[inds]==tab_net$pos1[ix_outlet]) ) ||
            any( (tab_net$pos2[inds]==tab_net$pos1[ix_outlet]) )
  linked_2= any( (tab_net$pos1[inds]==tab_net$pos2[ix_outlet]) ) ||
            any( (tab_net$pos2[inds]==tab_net$pos2[ix_outlet]) )
  rm(inds)
  if (sum(linked_1, linked_2) == 0) {
    stop("System outlet is not linked to upstream element(s).")
  }
  if (sum(linked_1, linked_2) > 1) {
    stop(paste("System outlet must not be linked to other element(s) on both ends. ",
      "Note: This also means that you cannot have a junction at the downstream end ",
      "of the reach defining the system's outlet. In such a case, you need to ",
      "manually edit the shape file and (re)-run TOPOCATCH with the edited version.",
      sep=""))
  }
  # Swap coordinates, if necessary (pos1 = upstream, pos2 = downstream)
  if (linked_2) {
    tmp= tab_net$pos1[ix_outlet]
    tab_net$pos1[ix_outlet]= tab_net$pos2[ix_outlet]
    tab_net$pos2[ix_outlet]= tmp
  }

  #.............................................................................
  # For each line, find the downstream neighbor. If neccessary, swap
  # end positions (if digitized in false direction).
  #.............................................................................

  if (!silent) print("Identifying downstream neighbor elements...")

  # Initialize ID of downstream neighbors to NA
  tab_net= cbind(tab_net, id_down= NA)
  # Assign a special ID to the outlet (make sure it is not used already)
  ID_END_NODE= "END_OF_SYSTEM"
  if (any(tab_net$id == ID_END_NODE))
    stop(paste("Reserved ID '",ID_END_NODE,"' is already in use.",sep=""))
  tab_net$id_down[ix_outlet]= ID_END_NODE
  change= TRUE
#########
## OLD CODE WHICH HAS BEEN SUBSTITUTED BY THE VECTOR CODE BELOW
#  while (change) {
#    change= FALSE
#    for (i in which(is.na(tab_net$id_down))) {
#      for (k in which(!is.na(tab_net$id_down))) {
#        # Connection exist, direction OK (pos1 = upstream, pos2 = downstream)
#        if (tab_net$pos2[i] == tab_net$pos1[k]) {
#          tab_net$id_down[i]= tab_net$id[k]
#          change= TRUE
#        # Connection exists but endpoints of upstream line need to be swapped
#        } else if (tab_net$pos1[i] == tab_net$pos1[k]) {
#          tab_net$id_down[i]= tab_net$id[k]
#          change= TRUE
#          tmp= tab_net$pos1[i]
#          tab_net$pos1[i]= tab_net$pos2[i]
#          tab_net$pos2[i]= tmp
#        }
#      }
#    }
#  }
#########
  while (change) {
    change= FALSE
    # Indices of sections whose downstream neighbor (is | is not) known
    inds_done= which(!is.na(tab_net$id_down))
    inds_todo= which(is.na(tab_net$id_down))
    # Find links
    ix_properDir= match(tab_net$pos2[inds_todo], tab_net$pos1[inds_done])  # Digitized in proper (downstream) direction
    ix_falseDir= match(tab_net$pos1[inds_todo], tab_net$pos1[inds_done])   # Digitized in false direction
    # Register links with proper directions
    ix_tmp= which(!is.na(ix_properDir))
    if (length(ix_tmp) > 0) {
      tab_net$id_down[inds_todo[ix_tmp]]= tab_net$id[inds_done[ix_properDir[ix_tmp]]]
      change=TRUE
    }
    # Register links with improper direction --> correct direction
    ix_tmp= which(!is.na(ix_falseDir))
    if (length(ix_tmp) > 0) {
      tab_net$id_down[inds_todo[ix_tmp]]= tab_net$id[inds_done[ix_falseDir[ix_tmp]]]
      change=TRUE
      # Adjust direction
      tmp= tab_net$pos1[inds_todo[ix_tmp]]
      tab_net$pos1[inds_todo[ix_tmp]]= tab_net$pos2[inds_todo[ix_tmp]]
      tab_net$pos2[inds_todo[ix_tmp]]= tmp
    }
    rm(inds_done, inds_todo, ix_properDir, ix_falseDir, ix_tmp)
  }
#########

  unlinked= tab_net$id[is.na(tab_net$id_down)] 
  if (length(unlinked) > 0 ) {
    stop(paste(length(unlinked)," elements could not be linked. The IDs follow: ",
      "'",paste(unlinked,collapse="', '"),"'",sep=""))
  }

  #.............................................................................
  # Add artificial end node
  #.............................................................................
  tab_net= rbind(tab_net, data.frame(
    id= ID_END_NODE,
    class= NA,
    pos1= tab_net$pos2[which(tab_net$id_down == ID_END_NODE)],
    pos2= NA,
    id_down= paste("DOWNSTREAM_OF_",ID_END_NODE,sep="")
  ))

  #.............................................................................
  # Add a catchment to all elements having a catchment
  #.............................................................................

  if (!silent) print("Creating catchment objects...")

  template= tab_net[(tab_net$class %in% classes_with_catchment),]
  tab_cat= data.frame(id= paste(prefix_catch,template$id,sep=""))
  # Check for generated duplicate IDs
  dups= which(tab_cat$id %in% tab_net$id)
  if (length(dups) > 0)
    stop(paste("Duplicate IDs generated. Please don't use the following feature",
      " ID(s) in the input shape file: '",paste(tab_cat$id[dups],collapse="', '"),"'.",sep=""))
  tab_cat$class= classname_catch
  tab_cat$pos1= NA
  # Catchment runoff enters river system DOWNSTREAM of corresponding reach/element
  tab_cat$pos2= template$pos2
  tab_cat$id_down= template$id_down
  tab_net= rbind(tab_net, tab_cat)
  rm(tab_cat)

  #.............................................................................
  # For each element, determine the list of upstream neighbors
  #.............................................................................

  if (!silent) print("Identifying upstream elements...")

  upElems= lapply(X=tab_net$id, FUN=uplist, tab_net$id_down, tab_net$id)
  tab_net= cbind(tab_net,
    nUp= unlist(lapply(upElems,length)),
    upElemIDs= unlist(lapply(upElems,paste,collapse=SEPCHAR))
  )
  tab_net$upElemIDs= as.character(tab_net$upElemIDs)
  rm(upElems)

  #.............................................................................
  # Delete head reaches and head minireaches that do not have a catchment upstream.
  # Note: This applies to head reaches, because the runoff of the corresponding
  # catchment enters the river system at the downstream end of the reach (i.e.
  # the reaches practically do not receive any input). It also applies to elements
  # belonging to the minireach class if they represent head reaches (those elements
  # never have a catchment by definition).
  #.............................................................................

  if (!silent) print("Removing unnecessary head elements...")

  old= tab_net$upElemIDs
  repeat {
    tab_net= subset(tab_net, ((!(tab_net$class %in% c(classname_reach,classname_minireach))) | (tab_net$nUp > 0)))
    # Now, update the list of upstream neighbors
    tab_net$nUp= NULL
    tab_net$upElemIDs= NULL
    upElems= lapply(X=tab_net$id, FUN=uplist, tab_net$id_down, tab_net$id)
    tab_net= cbind(tab_net,
      nUp= unlist(lapply(upElems,length)),
      upElemIDs= unlist(lapply(upElems,paste,collapse=SEPCHAR))
    )
    tab_net$upElemIDs= as.character(tab_net$upElemIDs)
    rm(upElems)

    if (identical(tab_net$upElemIDs, old)) {
      break
    }
    old= tab_net$upElemIDs
  }

  #.............................................................................
  # Create nodes where an element has upstream neighbors
  #.............................................................................

  if (!silent) print("Inserting node objects...")

  # Nodes are essentially required upstream of elements having > 1 upstream neighbor.
  # However, it is save to insert a node even if there is only 1 upstream neighbor
  # (node of order 1). This has the advantage that one can rely on the existence of
  # a node object at ALL interfaces of between other objects. This guarantees, for
  # example, that gage objects can be inserted at ANY of the interfaces.
  # Note: In a typical watershed model without special objects nodes of order 1
  #       do not exist in practice.
  inds= which(tab_net$nUp > 0)
  if (length(inds) > 0) {
    # Create nodes and set properties
    tab_nod= data.frame(
      id=paste(prefix_node,tab_net$id[inds],sep=""),
      class= classname_node,
      pos1= tab_net$pos1[inds],
      pos2= NA,
      id_down= tab_net$id[inds],
      nUp= tab_net$nUp[inds],
      upElemIDs= tab_net$upElemIDs[inds]
    )
    # Correct upstream ID of the elements downstream of the nodes
    tab_net$nUp[inds]= 1
    tab_net$upElemIDs[inds]= as.character(tab_nod$id[1:nrow(tab_nod)])
    # Redirect upstream elements into the node
    inds= which(tab_net$id_down %in% tab_net$id[inds])
    tab_net$id_down[inds]= paste(prefix_node,tab_net$id_down[inds],sep="")
  }
  rm(inds)

  # Check that all junctions where dissolved
  if (any(tab_net$nUp > 1)) stop("Failed to insert nodes.")
  # Check for generated duplicate IDs
  dups= which(tab_nod$id %in% tab_net$id)
  if (length(dups) > 0)
    stop(paste("Duplicate IDs generated. Please don't use the following feature",
      " ID(s) in the input shape file: '",paste(tab_nod$id[dups],collapse="', '"),"'.",sep=""))
  # Merge tables
  tab_net= rbind(tab_net, tab_nod)
  rm(tab_nod)

  #.............................................................................
  # Remove artificial end node from the system
  #.............................................................................
  tab_net= subset(tab_net, tab_net$id != ID_END_NODE)

  #.............................................................................
  # Insert gage objects at user-defined positions
  #.............................................................................

  # Note: The positions must always coincide with the position of a node.
  #       The gage object is then inserted just downstream of the node.

  tab_gages$id= as.character(tab_gages$id)
  if (nrow(tab_gages) > 0) {
    # Treat gage coordinates like the coordinates in the shape file
    tab_gages$x= round(tab_gages$x / coord_tol)
    tab_gages$y= round(tab_gages$y / coord_tol)
    # Add additional attributes to gages table
    tab_gages$class=classname_gage
    tab_gages$pos1= paste(as.character(tab_gages$x), as.character(tab_gages$y), sep=COORDSEP)
    tab_gages$pos2= NA
    tab_gages$id_down=NA               # to be set later
    tab_gages$nUp=NA                   # to be set later
    tab_gages$upElemIDs=NA             # to be set later
    for (i in 1:nrow(tab_gages)) {
      # Find the corresponding node by position
      k= which((tab_net$class == classname_node) & (tab_net$pos1 == tab_gages$pos1[i]))
      if (length(k) != 1)
        stop(paste("Failed to insert gage object '",tab_gages$id[i],"'. Found ",length(k),
          " node objects at the specified position (x: ",tab_gages$x[i]*coord_tol,", ",
          "y: ",tab_gages$y[i]*coord_tol,"). Maybe the tolerance value for the comparison of",
          " coordinates needs adjustment?",sep=""))
      # Set undefined info in gages table
      tab_gages$id_down[i]= tab_net$id_down[k]
      tab_gages$nUp[i]= 1
      tab_gages$upElemIDs[i]= tab_net$id[k]
      # Update the downstream element of the node
      tab_net$id_down[k]= tab_gages$id[i]
      # Update the upstream info for the element downstream of the gage
      # Note: If the gage is the last element in the system, there is no downstream object
      if (tab_gages$id_down[i] != ID_END_NODE) {
        k2= which(tab_net$upElemIDs == tab_net$id[k])  # can only be a single upstream element 
        if (length(k2) != 1)
          stop("Failed to update object downstream of gage.")
        tab_net$upElemIDs[k2]= tab_gages$id[i]
        rm(k2)
      }
      rm(k)
    }
    # Remove undesired fields
    tab_gages$x=NULL
    tab_gages$y=NULL
    # Check for duplicate IDs
    dups= which(tab_gages$id %in% tab_net$id)
    if (length(dups) > 0)
      stop(paste("Duplicate IDs detected after inserting gage objects. Please modify the following",
        " gage ID(s): '",paste(tab_gages$id[dups],collapse="', '"),"'.",sep=""))
    # Merge tables
    tab_net= rbind(tab_net, tab_gages)
    rm(tab_gages)
  }
 
  #.............................................................................
  # Create subclass info (include number of inputs, i.e. node order for nodes) 
  #.............................................................................
  tab_net$class= as.character(tab_net$class)
  tab_net$subclass= tab_net$class
  inds= which(tab_net$class == classname_node)
  tab_net$subclass[inds]= paste(tab_net$class[inds],tab_net$nUp[inds],sep="_n")
  rm(inds)

  #.............................................................................
  # Create output tables (object declaration and object linkage)
  #.............................................................................
  
  # Object declaration
  if (!silent) print("Creating object declaration table...")
  tab_objDecl= data.frame(
    object= tab_net$id,
    objectGroup= tab_net$subclass
  )

  # Init object linkage table
  if (!silent) print("Creating object linkage table...")

####### THIS SECTION WAS REPLACED BY THE VECTOR CODE BELOW
#  # Init table to full length
#  tab_objLink= data.frame(
#    targetObject= rep(NA, sum(tab_net$nUp) * nrow(tab_var)),
#    targetVariable= NA, sourceObject= NA, sourceVariable= NA, forwardType= NA
#  )
#  # Fill table
#  iStart= 1
#  for (i in which(tab_net$nUp > 0)) {
#    upElemIDs= unlist(strsplit(tab_net$upElemIDs[i], SEPCHAR, fixed=T))
#    if (length(upElemIDs) > 1) {
#      tarVar= paste(rep(tab_var$target, length(upElemIDs)),
#        sort(rep(1:length(upElemIDs),nrow(tab_var))), sep="_")
#    } else {
#      tarVar= tab_var$target
#    }
#    numRows= length(upElemIDs) * nrow(tab_var)
#    rng= iStart:(iStart+numRows-1)
#   tab_objLink$targetObject[rng]= as.character(rep(tab_net$id[i], numRows))
#    tab_objLink$sourceObject[rng]= as.character(sort(rep(upElemIDs, nrow(tab_var))))
#    tab_objLink$targetVariable[rng]= as.character(tarVar)
#    tab_objLink$sourceVariable[rng]= as.character(rep(tab_var$source, length(upElemIDs)))
#    tab_objLink$forwardType[rng]= rep(TRUE, numRows)
#    iStart= iStart + numRows
#  }
#######

  ## START VECTOR CODE
  # Init table to full length
  tab_objLink= data.frame(
    targetObject= rep(NA, sum(tab_net$nUp) * nrow(tab_var)),
    targetVariable= NA, sourceObject= NA, sourceVariable= NA, forwardType= NA
  )
  # Filter for objects with upstream neighbors
  tmp1= tab_net[(tab_net$nUp > 0),c("id","nUp","upElemIDs")]
  # Vector of all upstream neighbors of all objects
  tmp2= as.character(unlist(strsplit(paste(tmp1$upElemIDs, collapse=SEPCHAR),
    SEPCHAR, fixed=T)))
  # Vector of the corresponding target object IDs
  tmp3= rep(tmp1$id,tmp1$nUp)
  # Vector of variable name extensions
  tmp4= paste("_",sequence(tmp1$nUp),sep="")
  tmp4[rep((tmp1$nUp == 1),tmp1$nUp)]= ""

  # Populate table separately for each exchanged variable but sorted
  # by the ID of the target objects
  for (i in 1:nrow(tab_var)) {
    rng= seq(from=i, to=length(tmp2)*nrow(tab_var)-nrow(tab_var)+i, by=nrow(tab_var))
    tab_objLink$targetObject[rng]= as.character(tmp3)
    tab_objLink$sourceObject[rng]= as.character(tmp2)
    tab_objLink$targetVariable[rng]= paste(tab_var$target[i],tmp4,sep="")
    tab_objLink$sourceVariable[rng]= paste(tab_var$source[i])
    tab_objLink$forwardType[rng]= rep(TRUE, length(rng))
  }
  if (any(is.na(tab_objLink)))
    stop("Failed to create object linkage table.")
  ## END VECTOR CODE
  

  #.............................................................................
  # Result table with downstream neighbors
  #.............................................................................
  tab_downObj= tab_net[,c("id","id_down")]

  #.............................................................................
  # Return result tables in a list
  #.............................................................................

  return(list(
    objDecl= tab_objDecl,
    objLink= tab_objLink, 
    downObj= tab_downObj
  ))
}

