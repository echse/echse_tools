
################################################################################
# PART 1: Private utility functions first
################################################################################

#-------------------------------------------------------------------------------
# Helper functions for fitting of power models
sse= function(obs, mod) { sum((mod - obs)^2) } 
objFun= function(param, model, obs, predictor) { sse(obs, model(param, predictor)) } 
fit= function(paramInit, model, obs, predictor) {
  op= optim(par= paramInit, fn=objFun, gr=NULL, model, obs, predictor, method= c("Nelder-Mead"), control=list(maxit=1000))
  if (op$convergence != 0) {
    tf= tempfile()
    write.table(x=data.frame(obs=obs, predictor=predictor),file=tf,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
    stop(paste("No convergence (error code: ",op$convergence,", message: ",op$message,
      "). Initial parameters: [",paste(paramInit,collapse=","),
      "]. A snapshot of the data was saved in file '",tf,"'.",sep=""))
  }
  return(op$par)
}
powerModel = function(param, predictor) { param[["a"]] * predictor^param[["b"]] }
powerFit= function(paramInit, obs, predictor) { fit(paramInit, powerModel, obs, predictor) }

#-------------------------------------------------------------------------------
# Determines the cross-sections used in the interpolation and their weights
findNearXS= function(indexRCH, tabRCH, tabXS, miniFun) {
  # Check reach index
  if (!(indexRCH %in% 1:nrow(tabRCH)))
    stop("Invalid reach index.")
  # Check table columns
  if (!all(c("object","x1","y1","x2","y2","upstr_area") %in% names(tabRCH)))
    stop("Missing info in tabRCH.")
  if (!all(c("id","xCenter","yCenter","upstr_area") %in% names(tabXS)))
    stop("Missing info in tabXS.")
  # Check function
  if (!is.function(miniFun))
    stop("Argument 'miniFun' is not a registered function.")
  if (!identical(names(formals(miniFun)), c("dist","diffArea")))
    stop("Formal arguments of 'miniFun' must be 'dist' and 'diffArea'.")
  # Determine distance between target reach and all x-sections
  reachX= mean(tabRCH$x1[indexRCH],tabRCH$x2[indexRCH])
  reachY= mean(tabRCH$y1[indexRCH],tabRCH$y2[indexRCH])
  dists= sqrt((reachX-tabXS$xCenter)^2 + (reachY-tabXS$yCenter)^2)
  # Determine difference in upstream areas between target reach and all x-sections
  diffArea= abs(tabXS$upstr_area - tabRCH$upstr_area[indexRCH])
  # Find the nearest x-section whose chatchment is just smaller than that of the target reach
  ix1= which(tabXS$upstr_area <= tabRCH$upstr_area[indexRCH])
  if (length(ix1) == 0)
    stop(paste("Problem at reach '",tabRCH$object[indexRCH],"'. Upstream area is ",
      tabRCH$upstr_area[indexRCH],". No cross-section with smaller upstream area available.",sep=""))
  ix2= which.min(miniFun(dist=dists[ix1], diffArea=diffArea[ix1]))
  pos_small= ix1[ix2]
  # Find the nearest x-section whose chatchment is just larger than that of the target reach
  ix1= which(tabXS$upstr_area >= tabRCH$upstr_area[indexRCH])
  if (length(ix1) == 0)
    stop(paste("Problem at reach '",tabRCH$object[indexRCH],"'. Upstream area is ",
      tabRCH$upstr_area[indexRCH],". No cross-section with larger upstream area available.",sep=""))
  ix2= which.min(miniFun(dist=dists[ix1], diffArea=diffArea[ix1]))
  pos_large= ix1[ix2]
  # Return IDs and weights of the two cross-sections
  diff_small= abs(tabXS$upstr_area[pos_small] - tabRCH$upstr_area[indexRCH])
  diff_large= abs(tabXS$upstr_area[pos_large] - tabRCH$upstr_area[indexRCH])
  if (diff_small == 0) {
    weight_small= 1
    weight_large= 0
  } else if (diff_large == 0) {
    weight_small= 0
    weight_large= 1
  } else {
    weight_small= 1 - diff_small / (diff_small + diff_large)
    weight_large= 1 - diff_large / (diff_small + diff_large)
  }
  return(data.frame(stringsAsFactors=FALSE,
    id=tabXS$id[c(pos_small,pos_large)],
    weight=c(weight_small, weight_large),
    dist= dists[c(pos_small,pos_large)],                                     # only for docu purpose
    upstr_area= c(tabXS$upstr_area[pos_small], tabXS$upstr_area[pos_large])  # only for docu purpose
  ))
}

#-------------------------------------------------------------------------------
# Computes coefficients of a power model which is a weighted average of
# two power models fitted to two independend data sets
avgPowerCoeff= function(x1, y1, weight1, x2, y2, weight2, almostOne) {
  if (abs(sum(c(weight1,weight2)) - 1) > 1.e-06)
    stop(paste("Bad sum of weights. Weight_1: ",weight1," weight_2: ",weight2,".",sep=""))
  slopeEst= function(x,y) { (max(y)-min(y)) / (max(x)-min(x)) }
  if (weight1 >= almostOne) {
    return(powerFit(paramInit=c(a=slopeEst(x=x1,y=y1), b=1), obs=y1, predictor=x1))
  } else if (weight2 >= almostOne) {
    return(powerFit(paramInit=c(a=slopeEst(x=x2,y=y2), b=1), obs=y2, predictor=x2))
  } else {
    cof1= powerFit(paramInit=c(a=slopeEst(x=x1,y=y1), b=1), obs=y1, predictor=x1)
    cof2= powerFit(paramInit=c(a=slopeEst(x=x2,y=y2), b=1), obs=y2, predictor=x2)
    x= c(x1,x2)
    y= powerModel(cof1,x) * weight1 + powerModel(cof2,x) * weight2
    return(powerFit(paramInit=c(a=slopeEst(x=x,y=y), b=1), obs=y, predictor=x))
  }
}

#-------------------------------------------------------------------------------
# Computes reasonable 'flows of interest' based on a given maximum flow rate
flowRange= function(qmax) {
  fraction= 10000
  mini= qmax/fraction
  if (mini < 1) {
    llim= -1*ceiling(abs(log10(mini)))
  } else {
    llim= floor(abs(log10(mini)))
  }
  maxi= qmax
  if (maxi < 1) {
    ulim= -1*floor(abs(log10(maxi)))
  } else {
    ulim= ceiling(abs(log10(maxi)))
  }
  # Raw result vector
  tmp= 10^seq(from=llim, to=ulim, by=1)
  v= sort(c(tmp, tmp*2, tmp*5))
  # Cut off unnecessary elements
  inds_ex= which(v >= maxi)
  if (length(inds_ex) > 1)
    v= v[1:(length(v)-length(inds_ex)+1)]
  inds_ex= which(v <= mini)
  if (length(inds_ex) > 1)
    v= v[(1+length(inds_ex)-1):length(v)]
  return(v)
}

#-------------------------------------------------------------------------------
# Creates a data frame with the hydraulic characteristics of a reach.
reachHydrTable= function(
  roughness,  # see xs.steadyUniformWSE
  american,   # see xs.steadyUniformWSE
  slope,      # see xs.steadyUniformWSE
  fun_L2R,    # see xs.steadyUniformWSE
  fun_L2A,    # see xs.steadyUniformWSE
  par_L2R,    # see xs.steadyUniformWSE
  par_L2A,    # see xs.steadyUniformWSE
  reachlen,   # Length of the reach.
  q_max       # Highest flow rate to be tabulated. Lower bound = automatic.
) {
  # Check args
  checkArg(arg=roughness, len=1, type="numeric")
  checkArg(arg=american, len=1, type="logical")
  checkArg(arg=slope, len=1, type="numeric")
  if (slope <= 0)
    stop("Bed slope must be > 0.")
  checkArg(arg=fun_L2R, len=1, type="function")
  checkArg(arg=fun_L2A, len=1, type="function")
  checkArg(arg=reachlen, len=1, type="numeric")
  checkArg(arg=q_max, len=1, type="numeric")
  # Init output
  out= data.frame(Q= flowRange(q_max), D= NA, A= NA, V= NA, dVdQ= NA)
  # Set limits
  minWSE= 1.0e-9
  maxWSE= 200.
  # Loop over flow rates
  for (k in 1:nrow(out)) { 
    # Compute the normal depth for the current flow rate
    tryCatch({
      out$D[k]= xs.steadyUniformWSE(
        roughness=roughness, american=american, slope=slope,
        fun_L2R=fun_L2R, fun_L2A=fun_L2A, par_L2R=par_L2R, par_L2A=par_L2A,
        flow= out$Q[k], minWSE=minWSE, maxWSE=maxWSE, maxIter= 1000)
    }, error= function(e) {
      qAtMinWSE= xs.steadyUniformFlow(
        roughness=roughness, american=american, slope=slope, WSE=minWSE,
        fun_L2R=fun_L2R, fun_L2A=fun_L2A, par_L2R=par_L2R, par_L2A=par_L2A)
      qAtMaxWSE= xs.steadyUniformFlow(
        roughness=roughness, american=american, slope=slope, WSE=maxWSE,
        fun_L2R=fun_L2R, fun_L2A=fun_L2A, par_L2R=par_L2R, par_L2A=par_L2A)
      stop(paste("Failed to compute normal depth for Q=",out$Q[k],
        ". Tested boundary values: Q(WSE=",minWSE,")=",qAtMinWSE,", Q(WSE=",
        maxWSE,")=",qAtMaxWSE,".",sep=""))
    }) # End try
    # Compute wet area & storage volume at normal depth
    out$A[k]= fun_L2A(par_L2A, out$D[k])
    out$V[k]= out$A[k] * reachlen
  }
  # Compute derivative of storage with respect to flow
  # Part 1: At margins of the table
  n= nrow(out)
  out$dVdQ[1]= (out$V[2]-out$V[1]) / (out$Q[2]-out$Q[1])
  out$dVdQ[n]= (out$V[n]-out$V[n-1]) / (out$Q[n]-out$Q[n-1])
  # Part 2: Central rows of the table
  if (n > 3) {
    for (k in 2:(n-1)) {
      q_lower= (out$Q[k]+out$Q[k-1]) / 2
      dvdq_lower= (out$V[k]-out$V[k-1]) / (out$Q[k]-out$Q[k-1])
      q_upper= (out$Q[k+1]+out$Q[k]) / 2
      dvdq_upper= (out$V[k+1]-out$V[k]) / (out$Q[k+1]-out$Q[k])
      out$dVdQ[k]= approx(x=c(q_lower,q_upper), y=c(dvdq_lower,dvdq_upper),
        xout=out$Q[k])$y
    }
  }
  # Add record for zero flow
  out_q0= data.frame(Q= 0, D= 0, A= 0, V= 0, dVdQ= out$dVdQ[1])
  out= rbind(out_q0, out)
  # Round values in table
  out$D= round(out$D,3)
  out$A= round(out$A,3)
  out$V= round(out$V,3)
  out$dVdQ= round(out$dVdQ,2)
  # Return table
  return(out)
}

################################################################################
# PART 2: Exported methods
################################################################################

#' Estimation of reach properties for flood routing
#'
#' This function computes hydraulic characteristics for river reaches to be
#' used in flood-routing algorithms of hydrological models. In particular, the
#' output is useful for flood-routing techniques based on the approach of a
#' piece-wise linear reservoir. To use this function, cross-section data
#' must be available at least for some sites in a river basin. These sites are
#' not necessarily identical with the reaches of interest. See below for details.
#'
#' @param tabRCH A data frame listing basic properties of the reaches of interest.
#'   There must be at least the folling columns: 'object' (ID of reach object), 'x1', 'x2', 'y1',
#'   'y2' (coordinates of the reach's end-points), 'upstr_area'
#'   (size of upstream catchment), 'slope' (slope of river bed), 'roughness'
#'   (Manning's n or Strickler's k), 'length' (length of the reach).
#'   Note that a suitable data frame is returned in the
#'   \code{fileAttrRCH} argument of the \code{\link{hydroModelData}} method. The
#'   estimates of the slope in this data frame might need a plausibility check,
#'   however, as there is no guarantee that all values are greater than zero.
#' @param tabXS A data frame with the two columns 'datafile' and 'upstr_area'.
#'   Entries in the 1st column should point to files containing the geometry of
#'   individual cross-sections. The 2nd field must contain the catchment sizes
#'   corresponding to the cross-sections. These values are most efficiently
#'   identified in a GIS based on the shape file returned by the
#'   \code{\link{hydroModelData}} method.
#'   The data files referenced in the 1st column must be TAB-separated text files
#'   with (at least) the columns 'x', 'y', 'offset', and 'elevation'. Values in the offset
#'   column represent distances from a point at the bank (usually left bank) and
#'   the values must increase by row. Note that suitable data files are
#'   generated by the \code{\link{xs.extractDEM}} method, for example.
#' @param dirXS Name of a directory where the cross-section geometry files
#'   listed in \code{tabXS} reside. This is for the usual case that file
#'   names listed in \code{tabXS} are relative or just basenames.
#' @param dirOUT Name of a directory where output files are to be written to.
#' @param american A logical value. See the corresponding argument of the
#'   \code{\link{xs.steadyUniformFlow}} method for details.
#' @param fun_qMax A function which takes a catchment area as argument and
#'   returns the probable maximum flow rate. Note that the catchment is expected
#'   to in the same units as in the respective fields of \code{tabRCH} and
#'   \code{tabXS}. The function is used to restrict hydraulic computations to
#'   a reasonable range of flow rates. The default function expects the catchment
#'   size in units of square km. It is based on a global study by
#'   Herschy, R. (2001): The world's maximum observed floods, IAHS Publ., vol. 271. 
#' @param fun_hMax Similar to \code{fun_qMax} but this function must return a
#'   probable maximum flow depth instead of a flow rate. The default yields
#'   values of 2, 5, and 20 meters for catchment sizes of 1, 100, and 100.000
#'   square kilometers, respectively. This function is used to restrict the
#'   evaluation of the cross-section geometry data the approximate range of
#'   interest before fitting analytical functions to the data.
#' @param miniFun A function with two formal arguments 'dist' and 'diffArea'
#'   to control the selection of parent cross-sections for a particular reach
#'   from the available data pool.
#'   The 1st argument represents the spatial distance between the
#'   target reach and a potential parent cross-section. The 2nd argument is the
#'   corresponding absolute difference in the upstream catchment areas of the
#'   target reach and the potential parent cross-section, respectively.
#'   The parent cross-sections are selected so as to minimize the value of
#'   \code{miniFun}. The factor of 1000 in the default formulation accounts for
#'   the fact that the distance is in meters (if the x and y coordinates in the
#'   input tables are in meters) while the catchment areas are typically given
#'   in units of square kilometers.
#' @param almostOne A number > 0.5 and <= 1. If the interpolation weight for one
#'   of the two 'parent cross-sections' (see below) is greater than this value, this weight
#'   is assumed as 1 and the weight for the other cross-section is set to zero.
#'   Chosing a small value like 0.5 effectively turns the interpolation into a
#'   nearest neighbor approach.
#' @param replace If \code{TRUE}, existing files on \code{dirOUT} will be silently replaced.
#' @param plot If \code{TRUE}, a graphics files in PDF format is created for each
#'   target location. It can be used to visually assess the quality of the estimate.
#' @param silent If \code{TRUE}, some diagnostic messages are printed.
#'
#' @return The function's return value is \code{NULL}. The actual output is
#'   written to (usually a collection of) text file(s) containing
#'   all required in information for routing computations in the form of
#'   look-up tables. These tables are fully self-documenting.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de} while being at
#'   IIT Kharagpur, India
#'
#' @note To estimate the hydraulic properties of a river reach, its cross-section
#'   geometry must be know. In river basins, however, such cross-section are
#'   usually available for selected sites only. Therefore, in the first step,
#'   the functions 'interpolates' a cross-section for a reach, based on the data
#'   from two nearby 'parent cross-sections'. These two cross-sections are selected
#'   from the available data pool so that (1) the cross-sections are nearest to the reach of
#'   interest and (2) the 1st cross-section has a smaller and the 2nd one has a
#'   larger upstream catchment compared to the upstream catchment of the reach.
#'   A cross-section for the reach of interest is then estimated by linearly interpolating
#'   between the two parent cross-sections using the difference in the upstream catchment size
#'   as weights.
#'   
#'   This interpolation is, however, \emph{not} applied to the plain
#'   cross-section's geometry data but to derived characteristics, namely the
#'   relations \eqn{a(d)} and \eqn{r(d)} where \eqn{d} is the flow depth and
#'   \eqn{a} and \eqn{r} represent the wet area and the hydraulic radius,
#'   respectively. First, power models are fitted to represent those functions
#'   at the parent cross-sections. Second, 'intermediate' power models are
#'   determined by the above-mentioned interpolation to obtain estimates of
#'   \eqn{a(d)} and \eqn{r(d)} for the reach of interest. These two functions
#'   can then be used to compute the reach's hydraulic properties (including a
#'   rating curve, for example).
#'   
#'   The rationales behind this approach are: (1) Geology may be heteorogenious
#'   in larger river basins. Hence, it makes sense to assume that two nearby
#'   cross-sections are likely to be more 'similar' in terms of its properties
#'   than two cross-sections being located far away from each other. This is
#'   why the parent cross-sections are selected considering a minimum-distance
#'   criterion (see argument \code{miniFun}).
#'   (2) The cross-section geometry (namely the flow capacity) is assumed to be
#'   corellated with the size of the upstream catchment. Therefore, the
#'   parent cross-section whose catchment size is closer to the catchment size
#'   of the reach of interest is assigned a greater weight in the interpolation.
#'
#' @export

xs.reachPars= function(
  tabRCH,
  tabXS,
  dirXS,
  dirOUT,
  american,
  fun_qMax= function(sqkm) { 500 * sqkm^0.43 },
  fun_hMax= function(sqkm) { 2 * sqkm^0.2 },
  miniFun= function(dist, diffArea) { dist/1000 * diffArea },
  almostOne= 0.95,
  replace= FALSE,
  plot= FALSE,
  silent= TRUE
) {

  # Check args
  if (!silent) "Checking arguments..."
  checkArg(arg=tabRCH, len=NULL, type="data.frame")
  checkArg(arg=tabXS, len=NULL, type="data.frame")
  checkDir(dirXS, mustBeEmpty=FALSE)
  checkDir(dirOUT, mustBeEmpty=FALSE)
  checkArg(arg=american, len=1, type="logical")
  checkArg(arg=fun_qMax, len=1, type="function")  
  checkArg(arg=fun_hMax, len=1, type="function")
  checkArg(arg=miniFun, len=1, type="function")
  checkArg(arg=almostOne, len=1, type="numeric")
  if ((almostOne < 0.5) || (almostOne > 1))
    stop("Bad value for 'almostOne' argument.")
  checkArg(arg=replace, len=1, type="logical")
  checkArg(arg=silent, len=1, type="logical")

  # Check table with properties of target reaches
  if (!silent) print("Checking 1st input table...")
  cols= c("object","x1","x2","y1","y2","upstr_area","slope","roughness","length")
  miss= which(is.na(match(cols, names(tabRCH))))
  if (length(miss) > 0)
    stop(paste("Missing column(s) '",paste(cols[miss],collapse="', '"),"' in data frame 'tabRCH'.",sep=""))
  if (nrow(tabRCH) == 0)
    stop("Data frame 'tabRCH' is empty.")
  tabRCH$object= as.character(tabRCH$object)
  if (length(unique(tabRCH$object)) != nrow(tabRCH))
    stop("Entries in 'object' column of 'tabRCH' must be unique.")

  # Check table listing the catchment sizes of available cross-sections
  if (!silent) print("Checking 2nd input table...")
  cols= c("datafile","upstr_area")
  miss= which(is.na(match(cols, names(tabXS))))
  if (length(miss) > 0)
    stop(paste("Missing column(s) '",paste(cols[miss],collapse="', '"),"' in data frame 'tabXS'.",sep=""))
  if (nrow(tabXS) == 0)
    stop("Data frame 'tabXS' is empty.")
  if (length(unique(tabXS$datafile)) != nrow(tabXS))
    stop("File names listed in 'tabXS' must be unique.")

  # Read geometry data for all x-sections into a list of data frames
  # We also determine average coordinates for all x-sections to speed up nearest-neighbor search
  if (!silent) print("Reading geometry data...")
  xs= vector("list", length=nrow(tabXS))
  tabXS= cbind(tabXS, id=tabXS$datafile, xCenter=NA, yCenter=NA)
  for (i in 1:nrow(tabXS)) {
    f= paste(dirXS,"/",tabXS$datafile[i],sep="")
    if (!file.exists(f))
      stop(paste("Cross-section data file '",f,"' not found.",sep=""))
    geo= read.table(file=f, sep="\t", header=TRUE)
    cols= c("x","y","offset","elevation")
    miss= which(is.na(match(cols, names(geo))))
    if (length(miss) > 0)
      stop(paste("Missing column(s) '",paste(cols[miss],collapse="', '"),"' in cross-section data file '",f,"'.",sep=""))
    # Save this cross-section's geometry
    xs[[i]]= geo[,c("offset","elevation")]
    # Append average coordinates to table
    tabXS$xCenter[i]= mean(geo$x)
    tabXS$yCenter[i]= mean(geo$y)
    rm(geo)
  }

  if (!silent) print("Processing target reaches...")
  # Loop through target reaches
  for (indexRCH in 1:nrow(tabRCH)) {
   
    # Find the two cross-sections used in the interpolation and their weights
    sourceXS= findNearXS(indexRCH=indexRCH, tabRCH=tabRCH, tabXS=tabXS, miniFun=miniFun)
    if (nrow(sourceXS) != 2)
      stop(paste("Failed to select cross-sections for reach '",tabRCH$object[indexRCH],"'.",sep=""))

    # Compute the hydraulic properties of the two cross-sections
    tryCatch({
      index1= match(sourceXS$id[1], tabXS$id)  # Note: Elements in xs are in the same order as rows in tabXS
      offsets1= xs[[index1]]$offset
      elevations1= xs[[index1]]$elevation - min(xs[[index1]]$elevation)
      hydProp1= xs.geoFuncs(offsets=offsets1, elevations=elevations1,
        WSE=seq(from=0, to=fun_hMax(tabRCH$upstr_area[indexRCH]), length.out=100))
    }, error= function(e) {
      stop(paste("Failed to compute hydraulic properties of parent xs no. 1 (ID='",
        sourceXS$id[1],"') for reach '",tabRCH$object[indexRCH],"'.",sep=""))
    })
    tryCatch({
      index2= match(sourceXS$id[2], tabXS$id)  # Note: Elements in xs are in the same order as rows in tabXS
      offsets2= xs[[index2]]$offset
      elevations2= xs[[index2]]$elevation - min(xs[[index2]]$elevation)
      hydProp2= xs.geoFuncs(offsets=offsets2, elevations=elevations2,
        WSE=seq(from=0, to=fun_hMax(tabRCH$upstr_area[indexRCH]), length.out=100))
    }, error= function(e) {
      stop(paste("Failed to compute hydraulic properties of parent xs no. 2 (ID='",
        sourceXS$id[2],"') for reach '",tabRCH$object[indexRCH],"'.",sep=""))
    })

    # Fit weighted-average power models for wet area and rhyd
    tryCatch({
      cof_rhyd= avgPowerCoeff(
        x1=hydProp1$WSE, y1=hydProp1$rhyd, weight1=sourceXS$weight[1],
        x2=hydProp2$WSE, y2=hydProp2$rhyd, weight2=sourceXS$weight[2],
        almostOne=almostOne)
    }, error= function(e) {
      stop(paste("Failed to fit power model (hydraulic radius) for reach '",
        tabRCH$object[indexRCH],"' using parent xs '",sourceXS$id[1],"' and '",
        sourceXS$id[2],"'. Details: ",e,sep=""))
    })
    tryCatch({
      cof_area= avgPowerCoeff(
        x1=hydProp1$WSE, y1=hydProp1$area, weight1=sourceXS$weight[1],
        x2=hydProp2$WSE, y2=hydProp2$area, weight2=sourceXS$weight[2],
        almostOne=almostOne)
    }, error= function(e) {
      stop(paste("Failed to fit power model (wet area) for reach '",
        tabRCH$object[indexRCH],"' using parent xs '",sourceXS$id[1],"' and '",
        sourceXS$id[2],"'. Details: ",e,sep=""))
    })

    # Compute the hydraulic properties of the target reach
    out= reachHydrTable(
      roughness= tabRCH$roughness[indexRCH],
      american= american,
      slope= tabRCH$slope[indexRCH],
      fun_L2R= powerModel,
      fun_L2A= powerModel,
      par_L2R= cof_rhyd,
      par_L2A= cof_area,
      reachlen= tabRCH$length[indexRCH],
      q_max= fun_qMax(tabRCH$upstr_area[indexRCH])
    )
    # Output table for the target reach
    ofile= paste(dirOUT,"/",tabRCH$object[indexRCH],".tab",sep="")
    if (file.exists(ofile) && (!replace))
      stop(paste("File '",ofile,"' already exists.",sep=""))    
    head= paste("# Hydraulic properties of a reach for STEADY UNIFORM flow","\n",
      "#   Object ID:      '",tabRCH$object[indexRCH],"'","\n",
      "#   Upstream area: ",tabRCH$upstr_area[indexRCH],"\n",
      "#   Reach length:  ",tabRCH$length[indexRCH],"\n",
      "#   Bottom slope:  ",tabRCH$slope[indexRCH],"\n",
      "#   Roughness:     ",tabRCH$roughness[indexRCH],"\n",
      "# 1st / 2nd parent cross-section:","\n",
      "#   IDs:       '",sourceXS$id[1],"' / '",sourceXS$id[2],"'","\n",
      "#   weights:   ",signif(sourceXS$weight[1],3)," / ",signif(sourceXS$weight[2],3),"\n",
      "#   dist.s :   ",round(sourceXS$dist[1],1)," / ",round(sourceXS$dist[2],1),"\n",
      "#   up. areas: ",sourceXS$upstr_area[1]," / ",sourceXS$upstr_area[2],"\n",
      "# Descr. of columns:","\n",
      "#   Q:    Stream flow (L/T)","\n",
      "#   D:    Normal depth, i.e. max. flow depth (L)","\n",
      "#   A:    Wet x-section area (L^2)","\n",
      "#   V:    Storage volume (L^3)","\n",
      "#   dVdQ: Est. derivative dV/dQ (T)",
      sep="")
    write(x=head, file=ofile, sep="", ncolumns=length(head), append=FALSE)
    write(x= names(out), file=ofile, sep="\t", ncolumns=ncol(out), append=TRUE)
    write.table(x= out, file=ofile, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)

    # Generate plots
    if (plot) {
      pdf(file=paste(ofile,".pdf",sep=""), paper="a4r")
      layout(matrix(1:4,ncol=2,byrow=FALSE))
   
      # General info
      plot(0:1,0:1,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
      legend("left",bty="n",legend=c(
        paste("Target reach: '",tabRCH$object[indexRCH],"'",sep=""),
        paste("Parent xs1: '",sourceXS$id[1],"'",sep=""),
        paste("Parent xs2: '",sourceXS$id[2],"'",sep="")
      ))

      # Geometry of the parent cross sections
      range_offsets= range(offsets1, offsets2)
      range_elevations= range(elevations1, elevations2)
      plot(range_offsets, range_elevations, type="n", xlab="Offset", ylab="Elev. (normalized)")
      lines(offsets1, elevations1, col="blue")
      lines(offsets2, elevations2, col="red")
      legend("top", bty="n", lty=1, col=c("blue","red"),legend=c("xs1","xs2"))

      # Original hydraulic functions of parent cross-sections and interpolated
      # characteristics for the target reach
      range_wse= range(hydProp1$WSE, hydProp2$WSE)
      range_rhyd= range(hydProp1$rhyd, hydProp2$rhyd)
      range_area= range(hydProp1$area, hydProp2$area)
      wseSample=c(seq(from=0.1, to=10, by=0.1), 12, 15, 20, 25)

      plot(range_wse, range_rhyd, type="n", xlab="WSE", ylab="hyd. radius")
      points(hydProp1$WSE, hydProp1$rhyd, col="blue")
      points(hydProp2$WSE, hydProp2$rhyd, col="red")
      lines(wseSample, powerModel(cof_rhyd,wseSample), col="black")
      legend("topleft", bty="n", pch=c(1,1,NA), lwd=c(NA,NA,1), col=c("blue","red","black"),
        legend=c("xs1","xs2","Est."))

      plot(range_wse, range_area, type="n", xlab="WSE", ylab="wet area")
      points(hydProp1$WSE, hydProp1$area, col="blue")
      points(hydProp2$WSE, hydProp2$area, col="red")
      lines(wseSample, powerModel(cof_area,wseSample), col="black")
      legend("topleft", bty="n", pch=c(1,1,NA), lwd=c(NA,NA,1), col=c("blue","red","black"),
        legend=c("xs1","xs2","Est."))

      graphics.off()
    } # End of plotting

  } # End of loop over target reaches

  return(invisible(NULL))
}

