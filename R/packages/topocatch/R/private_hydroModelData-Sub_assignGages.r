
################################################################################

#' Identify all objects affecting a particular object
#'
#' This function solves the problem of finding all objects in a network that
#' affect a particular object of interest which is part of the network.
#'
#' @param tab_objDecl A data frame representing an object declaration table.
#'   Such a data frame must contain (at least) the two columns 'object' and 'objectGroup'.
#'   The \code{\link{hydroModelData}} method returns a suitable input (see its
#'   \code{fileObjDecl} argument).
#' @param tab_objLink A data frame representing an object linkage table.
#'   Such a data frame must contain (at least) the five columns
#'   'sourceObject', 'targetObject', 'sourceVariable', 'targetVariable', and
#'   'forwardType'. The \code{\link{hydroModelData}} method returns a suitable
#'    input (see its \code{fileObjLink} argument).
#' @param gageObjects A vector holding the names (IDs) of the objects of interest.
#'   This should be valid names in R (letters, digits, underscore).
#'   In hydrological applications, these are objects containing (or representing)
#'   stream gages. The supplied objects must exist in the 'object' column of
#'   \code{tab_objDecl} and in the 'targetObject' column of \code{tab_objLink}.

#' @return A data identical to \code{tab_objDecl} but with as many additional
#'   columns as there are elements in the vector \code{gageObjects}. These
#'   additional columns are of type integer. A value of 1 indicates, that the
#'   object of interest (specified in the column name) is affected by the
#'   object given in the 'object' column of the respective row. In models with
#'   only feed-forward interactions, this means that the object in the 'object'
#'   column is located upstream of the object of interest (given as column name).
#'   In models with feedbacks, this is not necessarily so.
#'   A value of 0 indicates that the object of interest (column name) is not
#'   affected by the object given in the 'object' column.
#'
#' @note The function is applicable to dendritic (i.e. tree-like) systems as
#'   well as to meshed systems (i.e. systems with flow diversions).
#'
#'  The function is also capable of handling both types of object links, i.e.
#'  forward (upstream to downstream) as well as backward interactions (feedbacks).
#'
#'  The current implementation may be slow if applied to very large systems.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Generation of sample object declaration table (with spatial coordinates)
#' decl= data.frame(stringsAsFactors=FALSE,
#'   object=c("a","b","c","d","e","f","g","h","i","a2","b2","p","q","r","s","t","p2"),
#'   objectGroup="-",
#'   xpos= c(  1,  1,  2,  2,  1,  1,  2,  2,  1,   3,   3,  6,  5,  4,  3,  3,   6),
#'   ypos= c(  9,  8,  7,  6,  5,  4,  3,  2,  1,   9,   8,  8,  7,  6,  5,  4,   6)
#' )
#' # Generation of corresponding object linkage table
#' link= data.frame(stringsAsFactors=FALSE,
#'   sourceObject=c("a","b","c","d","d","e","f","g","h","a2","b2","p","q","r","s","t","p2","t"),
#'   targetObject=c("b","c","d","e","s","f","g","h","i","b2","c", "q","r","s","t","g","q","d"),
#'   sourceVariable="-", targetVariable="-",
#'   forwardType=as.logical(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0))
#' )
#' # Definition of gage objects
#' gages= c("s", "h", "e", "r", "b")
#'
#' # Call method 
#' result=assignGages(tab_objDecl=decl[,1:2], tab_objLink=link, gageObjects=gages)
#' # Plot basic network
#' plot(range(decl$xpos)+c(-1,1), range(decl$ypos)+c(-1,1),
#'   type="n", xlab="x", ylab="y")
#' for (i in 1:nrow(link)) {
#'   i_src= match(link$sourceObject[i], decl$object)
#'   i_tar= match(link$targetObject[i], decl$object)
#'   lines(x=decl$xpos[c(i_src,i_tar)], y=decl$ypos[c(i_src,i_tar)], col="grey",
#'     lty=2-link$forwardType[i])
#' }
#' points(decl$xpos, decl$ypos, col="grey", cex=3)
#' text(decl$xpos, decl$ypos, decl$object)
#' # Show gages and their upstream objects
#' clr= rainbow(length(gages))
#' for (i in 1:length(gages)) {
#'   points(x=decl$xpos[decl$object==gages[i]], y=decl$ypos[decl$object==gages[i]],
#'     col=clr[i], pch=20, cex=3)
#'   text(x=decl$xpos[decl$object==gages[i]], y=decl$ypos[decl$object==gages[i]],
#'     gages[i])
#'   gagedObj= result$object[which(result[,gages[i]] != 0)]
#'   if (length(gagedObj) > 0) {
#'     sub= subset(decl, decl$object %in% gagedObj)
#'     points(sub$xpos, sub$ypos, cex=3+i, col=clr[i])
#'   }
#' }
#' legend("bottomright", bty="n", pch=c(20,1), legend=c("Gage object","Influencing obj."))
#' legend("topright", bty="n", lty=c(1,2), legend=c("Forward link","Backward link"))

assignGages= function(
  tab_objDecl,
  tab_objLink,
  gageObjects
) {
  # Check args
  checkArg(arg=tab_objDecl, len=NULL, type="data.frame")
  checkArg(arg=tab_objLink, len=NULL, type="data.frame")
  checkArg(arg=gageObjects, len=NULL, type="character")

  # Check tables
  colnames= c("object","objectGroup")
  bad= which(!(colnames %in% names(tab_objDecl)))
  if (length(bad) > 0)
    stop(paste("Missing column(s) '",paste(colnames[bad],collapse="', '"),
      "' in object declaration table.",sep=""))
  colnames= c("targetObject","targetVariable","sourceObject","sourceVariable","forwardType")
  bad= which(!(colnames %in% names(tab_objLink)))
  if (length(bad) > 0)
    stop(paste("Missing column(s) '",paste(colnames[bad],collapse="', '"),
      "' in object linkage table.",sep=""))

  # Check supplied gage objects
  if (length(gageObjects) == 0)
    stop("No gage objects specified.")
  # Values must be unique
  if (length(gageObjects) != length(unique(gageObjects)))
    stop("Duplicate gage objects detected.")
  # Objects must exist in both the tables
  bad= which(!(gageObjects %in% tab_objDecl$object))
  if (length(bad) > 0)
    stop(paste("Gage object(s) '",paste(gageObjects[bad],collapse="', '"),
      "' not existent in object declaration table.",sep=""))
  bad= which(!(gageObjects %in% tab_objLink$targetObject))
  if (length(bad) > 0)
    stop(paste("Gage object(s) '",paste(gageObjects[bad],collapse="', '"),
      "' not existent as target objects in object linkage table.",sep=""))

  # Simplify info in object linkage table
  link= tab_objLink
  rm(tab_objLink)
  link$targetVariable= NULL
  link$sourceVariable= NULL
  link$forwardType= NULL
  link= link[!duplicated(link),] # there may be multi-variable relations between two particular elements
  names(link)[names(link) == "sourceObject"]= "src"
  names(link)[names(link) == "targetObject"]= "tar"
  link$tar= as.character(link$tar)
  link$src= as.character(link$src)

  # Initialize matrix with gage info (nrow= number of links; ncol= number of gages)
  # --> Value FALSE means that the gage is not affected by this link
  gage= matrix(FALSE, ncol=length(gageObjects), nrow=nrow(link))
  tmp= match(link$tar, gageObjects)       # treat the target as the gage object
  for (i in 1:length(tmp)) {
    if (!is.na(tmp[i]))
      gage[i, tmp[i]]= TRUE
  }

  # Note: (1) One target object may have multiple source objects --> Junctions
  #       (2) One source object may have multiple target objects --> Flow diversions (split-flow junctions)
  # 
  # Therefore, we cannot simply walk along flow paths but must use an iterative approach that
  # treats all links separately.

  # For each link, identify the affected gages
  change= TRUE
  while (change) {
    change= FALSE
    for (i in 1:nrow(link)) {
      # Find all upstream neighbor links
      k= which(link$tar == link$src[i])
      if (length(k) > 0) {
        # Does upstream neighbor already have (at least) the info of the downstream neighbor?
        tmp= which(gage[i,])
        if (length(tmp) > 0) {
          if (!all(gage[k,tmp]))
            change=TRUE
          # Transfer info of upstream neighbor
          gage[k,tmp]= TRUE
        }
      }
    }
  }

  # Collapse the result for rows with the same source object. In this case, the
  # source elements are the origin of a diversion (flow split) and the info on
  # all downstream branches needs to be collected.
  effects= aggregate(x=gage, by=list(sourceObject=link$src), FUN=any)
  names(effects)[2:ncol(effects)]= gageObjects

  # Let the gage objects also affect themselves. Advantage: The gage of interest
  # is included in the list of objects when the object declaration table is filtered
  # for the catchment of this gage
  for (g in gageObjects) {
    effects[effects$sourceObject==g,g]= TRUE
  }

  # Convert info from logical to integer (smaller output file)
  for (i in 2:ncol(effects)) {
    effects[,i]= as.integer(effects[,i])
  }

  # Create output
  result= merge(x= tab_objDecl, y=effects, by.x="object", by.y="sourceObject", all=TRUE)
  if (nrow(result) != nrow(tab_objDecl))
    stop("Bad length of result table.")
  result[is.na(result)]= 0

  return(result)

}

