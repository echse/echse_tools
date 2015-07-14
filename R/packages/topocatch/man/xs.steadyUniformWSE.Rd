\name{xs.steadyUniformWSE}
\alias{xs.steadyUniformWSE}
\title{Water surface elevation under steady uniform flow}
\usage{
  xs.steadyUniformWSE(roughness, american, slope, fun_L2R,
    fun_L2A, par_L2R, par_L2A, flow, minWSE, maxWSE,
    maxIter = 1000)
}
\arguments{
  \item{flow}{The flow rate for which the water surface
  elevation is to be computed. Must be a scalar.}

  \item{minWSE}{Lower boundary of the search interval.
  Should be set to a water surface elevation where flow is
  known to be zero. Typically, this is an elevation equal
  to or lower than the lowest point in a cross-section.}

  \item{maxWSE}{Upper boundary of the search interval.
  Should be set to the highest thinkable water surface
  elevation for the cross-section. Such a value can
  typically be guessed considering the range of reasonable
  flow depths.}

  \item{maxIter}{The maximum allowed number of iterations.}

  \item{roughness}{The 'roughness' value to account for
  energy losses due to friction and turbulence. Must be a
  scalar. See \code{american} to possible values.}

  \item{american}{A logical value. If \code{TRUE},
  \code{roughness} is interpreted as Manning's \emph{n}
  with typical values being in range 0.01 to 0.2. If
  \code{FALSE}, the values of \code{roughness} should be in
  the approximate range 5 to 100 and represent Strickler's
  \emph{k} (\eqn{k=1/n}).}

  \item{slope}{The slope of the channel's bed as a
  dimensionless number (e.g. m/m).}

  \item{fun_L2R}{A function to return the cross-secion's
  hydraulic radius for the values in vector \code{WSE}. The
  interface should be \code{f(params,wse)}. When the
  function is called, the value of \code{par_L2R} is passed
  as the first argument and the vector \code{WSE} is passed
  as the second one.}

  \item{fun_L2A}{A function to return the cross-section's
  wer area for the values in vector \code{WSE}. The
  interface should be \code{f(params,wse)}. When the
  function is called, the value of \code{par_L2A} is passed
  as the first argument and the vector \code{WSE} is passed
  as the second one.}

  \item{par_L2R}{A second \emph{constant} argument (i.e.
  parameter) to be passed to \code{fun_L2R}. Use a list do
  pass more complex data items.}

  \item{par_L2A}{A second \emph{constant} argument (i.e.
  parameter) to be passed to \code{fun_L2A}. Use a list do
  pass more complex data items.}
}
\value{
  The water surface elevation corresponding to \code{flow}
  under steady, uniform conditions.
}
\description{
  The function computes the water surface elevation
  corresponding to a flow rate for steady, uniform
  conditions. This is the water surface elevation at the
  so-called 'normal depth'. The calculation is based
  Manning's equation (also known as
  Gaukler-Manning-Strickler's eqn.). The
  \code{\link{uniroot}} method is used to solve the
  root-finding problem.
}
\note{
  An error is generated if no solution is found within
  \code{maxIter} iterations. Note that, in contrast to the
  inverse function \code{\link{xs.steadyUniformFlow}}, this
  function expects a scalar argument, not a vector.
}
\examples{
# Generate cross-section data
geo= data.frame(x= 0:5, z= c(2,1.5,0.5,0.5,1.5,2))
# Wrapper functions for computing wet area and hydr. radius from data
L2A= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$area) }
L2R= function(geo, w) { return(xs.geoFuncs(geo$x, geo$z, w)$rhyd) }
print(xs.steadyUniformWSE(roughness=0.03, american=TRUE, slope=0.01,
  fun_L2R=L2R, fun_L2A=L2A, par_L2R=geo, par_L2A=geo, flow=0.3358317,
  minWSE=min(geo$z), maxWSE=max(geo$z), maxIter= 1000))
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\seealso{
  \code{\link{xs.steadyUniformFlow}} for the reverse
  computation
}

