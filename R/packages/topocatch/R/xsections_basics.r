
################################################################################
# Functions to compute hydraulic propierties of for river cross-sections
################################################################################

#' Basic hydraulic properties of a river cross-section
#'
#' The function returns values of the wet area, wet perimeter, and hydraulic
#' radius of a cross-section corresponding to a vector of given water
#' surface elevations (WSE).
#'
#' @param offsets A numeric vector of offsets. Typically, these are horizontal
#'   distances measured from a fixed point at one of the river banks. It does
#'   not matter whether the fix point is at the left or right bank.
#' @param elevations A numeric vector of elevations corresponding to the
#'   offsets in \code{offsets}. These can be absolute elevations (m a.s.l.) or
#'   vertical distances to any other convenient reference elevation.
#' @param WSE A vector of water surface elevations for which wet area and wet
#'   parameter are to be computed. These values are typically within
#'   \code{range{elevations}}.
#'
#' @return A list with the following components, each being a vector of
#'   \code{length(WSE)}.
#'   \item{WSE}{The input vector \code{WSE}.}
#'   \item{area}{Values of wet cross-section area.}
#'   \item{peri}{Values of wer perimeter.}
#'   \item{rhyd}{Values of hydraulic radius.}
#'
#' @note For values in \code{WSE} being smaller than \code{min(elevations)},
#'   all cross-section characteristics are zero. For values that exceede
#'   \code{max(elevations)}, the function returns valid results assuming
#'   vertical boundaries at the smallest and largest offset, respectively.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @export
#'
#' @examples
#' # Generate cross-section data
#' geo= data.frame(x= 0:5, z= c(2,1.5,0.5,0.5,1.5,2))
#' wsElevs= seq(from=0.5, to=1.75, by=0.25)
#' hydrProp= xs.geoFuncs(offsets=geo$x, elevations=geo$z, WSE=wsElevs)
#' plot(geo$x, geo$z, type="l", xlab="offset", ylab="elevation", asp=1)
#' abline(h=wsElevs, lty=3)
#' print(as.data.frame(hydrProp))

xs.geoFuncs= function(
  offsets,
  elevations,
  WSE=seq(from=min(elevations), to=max(elevations), length.out=5)
) {
  checkArg(arg=offsets, len=NULL, type="numeric")
  checkArg(arg=elevations, len=NULL, type="numeric")
  checkArg(arg=WSE, len=NULL, type="numeric")
  if (length(offsets) != length(elevations))
    stop("Offset and elevation vector are of different length.")
  if (length(offsets) < 2)
    stop("Not a valid cross-section. Too few offsets.")
  if (!identical(offsets, sort(offsets)))
    stop("Vector of offsets is not in increasing order.")
  if (length(WSE) < 1)
    stop("Number of elements in WSE must be 1 at least.")
  area= rep(0,length(WSE))
  peri= rep(0,length(WSE))
  rhyd= rep(0,length(WSE))
  for (i in 1:length(WSE)) {
    for (k in 1:(length(offsets)-1)) {
      minelev= min(elevations[k],elevations[k+1])
      maxelev= max(elevations[k],elevations[k+1])
      offsetdiff= abs(offsets[k]-offsets[k+1])
      # Bottom is dry
      if (WSE[i] <= minelev) {
        a= 0
        p= 0
      # Bottom completely submerged
      } else if (WSE[i] >= maxelev) {
        t= offsetdiff
        a= ((WSE[i] - elevations[k]) + (WSE[i] - elevations[k+1])) / 2 * t
        p= sqrt((elevations[k]-elevations[k+1])^2 + t^2)
      # Bottom partially submerged
      } else {
        t= offsetdiff * ((WSE[i] - minelev) / (maxelev - minelev))
        a= t * (WSE[i] - minelev) / 2
        p= sqrt((WSE[i]-minelev)^2 + t^2)
      }
      area[i]= area[i] + a
      peri[i]= peri[i] + p
    }
    if (peri[i] > 0) {
      rhyd[i]= area[i] / peri[i]
    } else {
      rhyd[i]= 0
    }
  }
  return(list(WSE= WSE, area=area, peri=peri, rhyd=rhyd))
}

################################################################################

#' Flow rate for steady uniform conditions
#'
#' The function computes flow rates corresponding to given water surface
#' elevations (WSE) under steady, uniform conditions. The calculation is based
#' Manning's equation (also known as Gaukler-Manning-Strickler's eqn.).
#'
#' @param roughness The 'roughness' value to account for energy losses due to
#'   friction and turbulence. Must be a scalar. See \code{american} to possible values.
#' @param american A logical value. If \code{TRUE}, \code{roughness}
#'   is interpreted as Manning's \emph{n} with typical values being in range
#'   0.01 to 0.2. If \code{FALSE}, the values of \code{roughness} should be in
#'   the approximate range 5 to 100 and represent Strickler's \emph{k} (\eqn{k=1/n}).
#' @param slope The slope of the channel's bed as a dimensionless number (e.g. m/m).
#' @param WSE A vector of water surface elevations for which the the flow rates
#'   are sought.
#' @param fun_L2R A function to return the cross-secion's hydraulic
#'   radius for the values in vector \code{WSE}. The interface should be
#'   \code{f(params,wse)}. When the function is called, the value of
#'   \code{par_L2R} is passed as the first argument and the vector \code{WSE} is
#'   passed as the second one.
#' @param fun_L2A A function to return the cross-section's wer area for
#'   the values in vector \code{WSE}. The interface should be
#'   \code{f(params,wse)}. When the function is called, the value of
#'   \code{par_L2A} is passed as the first argument and the vector \code{WSE} is
#'   passed as the second one.
#' @param par_L2R A second \emph{constant} argument (i.e. parameter) to be passed
#'   to \code{fun_L2R}. Use a list do pass more complex data items.
#' @param par_L2A A second \emph{constant} argument (i.e. parameter) to be passed
#'   to \code{fun_L2A}. Use a list do pass more complex data items.
#'
#' @return A vector of flow rates of the same length as \code{WSE} if
#'   all values are finite and positive. Otherwise, an error is generated and
#'   an appropriate message is issued.
#'
#' @note Use of the functions \code{fun_L2R} and \code{fun_L2A} to represent
#'   the cross-section makes this function more generic. It can be applied to
#'   cross-sections where the actual geometry data is available. It can also be
#'   applied, however, if the actual geometry is not available but the two
#'   characteristic functions are know (typically as power laws).
#'   See the examples for how to apply this function in the case where the
#'   actual geometry is known using calls to \code{\link{xs.geoFuncs}}.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @seealso \code{\link{xs.steadyUniformWSE}} for the reverse computation
#'
#' @export
#'
#' @examples
#' # Generate cross-section data and plot
#' geo= data.frame(x= 0:5, z= c(2,1.5,0.5,0.5,1.5,2))
#' wsElevs= seq(from=0.5, to=1.75, by=0.25)
#' hydrProp= xs.geoFuncs(offsets=geo$x, elevations=geo$z, WSE=wsElevs)
#' plot(geo$x, geo$z, type="l", xlab="offset", ylab="elevation", asp=1)
#' abline(h=wsElevs, lty=3)
#' # Wrapper functions for computing wet area and hydr. radius from data
#' L2A= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$area) }
#' L2R= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$rhyd) }
#' # Compute flow rates
#' hydrProp$flow= xs.steadyUniformFlow(roughness=0.03, american=TRUE, slope=0.01,
#'   WSE=wsElevs, fun_L2R=L2R, fun_L2A=L2A, par_L2R=geo, par_L2A=geo)
#' print(as.data.frame(hydrProp))
#' # Note: This example may performance sub-optimal for large amounts of
#' #       data because xs.geoFuncs is called twice in xs.steadyFlow.

xs.steadyUniformFlow= function(
  roughness,
  american,
  slope,
  WSE,
  fun_L2R,
  fun_L2A,
  par_L2R,
  par_L2A
) {
  # Check args
  checkArg(arg=roughness, len=1, type="numeric")
  checkArg(arg=american, len=1, type="logical")
  if (american && (roughness > 1))
    stop(paste("Roughness of '",roughness,"' not reasonable for american==TRUE.",sep=""))
  if ((!american) && (roughness < 1))
    stop(paste("Roughness of '",roughness,"' not reasonable for american==FALSE.",sep=""))
  checkArg(arg=slope, len=1, type="numeric")
  checkArg(arg=WSE, len=NULL, type="numeric")
  if (length(WSE) < 1)
    stop("Empty vector of water surface elevations.")
  checkArg(arg=fun_L2R, len=1, type="function")
  checkArg(arg=fun_L2A, len=1, type="function")
  # Compute
  if (american) roughness= 1/roughness
  result= roughness * sqrt(slope) * fun_L2A(par_L2A, WSE) * fun_L2R(par_L2R, WSE)^(2/3)
  # Report problems
  bad= which((!is.finite(result)) | (result < 0))
  if (any(bad)) {
    stop(paste("Manning's equation returns bad value(s) for the WS elevation(s) ",
    paste(WSE,collapse=","),". The corresponding values of wet area and ",
    "hydraulic radius are ",paste(fun_L2A(par_L2A, WSE),collapse=","),
    " and ",paste(fun_L2R(par_L2R, WSE),collapse=","),". Other used inputs were ",
    "Roughness: ",roughness," American?: ",american," Slope: ",slope,".",sep=""))
  }
  return(result)
}

################################################################################

# Helper function for estimation of normal depth --> Root of this function is sought
flowDiff= function(WSE, roughness, american, slope, fun_L2R, fun_L2A, par_L2R, par_L2A, flow) {
  return(
    xs.steadyUniformFlow(roughness, american, slope, WSE, fun_L2R, fun_L2A, par_L2R, par_L2A) - flow
  )
}

#' Water surface elevation under steady uniform flow
#'
#' The function computes the water surface elevation corresponding to a flow rate for
#' steady, uniform conditions. This is the water surface elevation at the so-called
#' 'normal depth'. The calculation is based Manning's equation
#' (also known as Gaukler-Manning-Strickler's eqn.). The \code{\link{uniroot}}
#' method is used to solve the root-finding problem.
#'
#' @inheritParams xs.steadyUniformFlow
#' @param flow The flow rate for which the water surface elevation is to be
#'   computed. Must be a scalar.
#' @param minWSE Lower boundary of the search interval. Should be set to a
#'   water surface elevation where flow is known to be zero. Typically, this is
#'   an elevation equal to or lower than the lowest point in a cross-section.
#' @param maxWSE Upper boundary of the search interval. Should be set to the
#'   highest thinkable water surface elevation for the cross-section. Such a
#'   value can typically be guessed considering the range of reasonable flow depths.
#' @param maxIter The maximum allowed number of iterations.
#'
#' @return The water surface elevation corresponding to \code{flow} under
#'   steady, uniform conditions.
#'
#' @author David Kneis \email{david.kneis@@uni-potsdam.de}
#'
#' @note An error is generated if no solution is found within \code{maxIter}
#'   iterations. Note that, in contrast to the inverse function
#'   \code{\link{xs.steadyUniformFlow}}, this function expects a scalar argument,
#'   not a vector.
#'
#' @seealso \code{\link{xs.steadyUniformFlow}} for the reverse computation
#'
#' @export
#'
#' @examples
#' # Generate cross-section data
#' geo= data.frame(x= 0:5, z= c(2,1.5,0.5,0.5,1.5,2))
#' # Wrapper functions for computing wet area and hydr. radius from data
#' L2A= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$area) }
#' L2R= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$rhyd) }
#' print(xs.steadyUniformWSE(roughness=0.03, american=TRUE, slope=0.01,
#'   fun_L2R=L2R, fun_L2A=L2A, par_L2R=geo, par_L2A=geo, flow=0.3358317,
#'   minWSE=min(geo$z), maxWSE=max(geo$z), maxIter= 1000))

xs.steadyUniformWSE= function(
  roughness,
  american,
  slope,
  fun_L2R,
  fun_L2A,
  par_L2R,
  par_L2A,
  flow,
  minWSE,
  maxWSE,
  maxIter= 1000
) {
  # Check args
  checkArg(arg=roughness, len=1, type="numeric")
  checkArg(arg=american, len=1, type="logical")
  checkArg(arg=slope, len=1, type="numeric")
  checkArg(arg=fun_L2R, len=1, type="function")
  checkArg(arg=fun_L2A, len=1, type="function")
  checkArg(arg=flow, len=NULL, type="numeric")
  if (length(flow) < 1)
    stop("Empty vector of flow rates.")
  checkArg(arg=minWSE, len=1, type="numeric")
  checkArg(arg=maxWSE, len=1, type="numeric")
  if (minWSE >=  maxWSE)
    stop("Bad value for minWSE or maxWSE.")
  checkArg(arg=maxIter, len=1, type="integer")
  # Solve
  tryCatch({
    res= uniroot(f=flowDiff, interval=c(minWSE, maxWSE),
      roughness=roughness, american=american, slope=slope,
      fun_L2R=fun_L2R, fun_L2A=fun_L2A, par_L2R=par_L2R, par_L2A=par_L2A, flow=flow,
       maxiter=maxIter)
  }, error= function(e) {
    stop(paste("Computation of water surface elevation failed. The search range for",
      " root-finding was ",minWSE," to ",maxWSE,". Error details: ",e,sep=""))
  })
  if (res$iter >= maxIter)
    stop("No convergence in computation of water surface elevation.")
  return(res$root)
}


