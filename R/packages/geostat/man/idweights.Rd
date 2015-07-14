\name{idweights}
\alias{idweights}
\title{Compute weights for inverse-distance spatial interpolation}
\usage{
  idweights(x_tar, y_tar, x_src, y_src, nsectors = 4,
    norigins = 1, power = 2)
}
\arguments{
  \item{power}{The power to be used when computing inverse
  distance weights.  In many applications, a value of 2 is
  used as a default. Higher values lead to increased
  weights for nearer points. A value of zero makes the
  weights independend of the points' distance and results
  in simple averaging.}

  \item{x_tar}{x-coordinate of the target point (scalar)}

  \item{y_tar}{y-coordinate of the target point (scalar)}

  \item{x_src}{x-coordinates of the source points (vector)}

  \item{y_src}{y-coordinates of the source points (vector)}

  \item{nsectors}{The number of sectors to use (minimum 1).
  Each sector will be 360 / \code{nsectors} degrees wide.
  Setting \code{nsectors} to 1 results in a
  nearest-neighbor approach (also known as Thiessen
  method).}

  \item{norigins}{The number of sector origins to be
  testet. If \code{norigins} is > 1, the search sectors
  will be rotated with an increment of 90 / \code{norigins}
  degrees around the center point given by \code{x_tar} and
  \code{y_tar}. The returned results correspond to the best
  angle of origin tested. The criterias of optimality are
  the number of non-empty sectors and the average distance
  between the target point and the source points. If
  \code{nsectors} == 1, sector rotation does not make sense
  and \code{norigins} = 1 should be used. In a standard
  setting with \code{nsectors} ==4, one would typically set
  \code{norigins} to a value in range 1 through 4. Higher
  values may give better results at the cost of increased
  computation time.}
}
\value{
  A numeric vector of the same length as \code{x_src}
  holding the inverse-distance interpolation weights. The
  order of elements is the same as in \code{x_src} and
  \code{y_src}.
}
\description{
  Given a point 'T' with coordinates \code{x_tar, y_tar}
  and a set of points 'S' defined by the coordinate vectors
  \code{x_src} and \code{y_src}, the function computes
  interpolation weights for all point in 'S'. These weights
  can be used to estimate the an unknown value of a
  variable at point 'T' from known values at the points
  'S'. The weights are computed by the inverse-distance
  approach, i.e. the value of the weight decreases as the
  distance between 'T' and the particular point from set
  'S' increases.
}
\note{
  The selection of suitable points from set 'S' is
  performed by the \code{\link{neighbors}} method. A weight
  of zero is assigned to all points from 'S' not being part
  of the selection (because preference was given to a
  nearer point in the same sector).

  Use of the return value: The returned vector can be used
  to filter the set of source points. For that purpose, the
  return vector needs to be converted into a logical mask
  by checking which elements differ significantly from zero
  using something like \code{mask = x > 0.01}.
}
\examples{
# Demonstrates the regionalization of point observations to a regular grid
nsources= 16
xrng= c(0,100)
yrng= c(0,100)
dx= 4
dy= 4
nsectors= 12
norigins= 10
power= 2.
ncolors= 100
# Generate some source data
sources= data.frame(
  x= min(xrng) + runif(nsources) * diff(xrng),
  y= min(yrng) + runif(nsources) * diff(yrng),
  value= rnorm(mean=0, sd=1, n=nsources)
)
# Generate target locations (regular grid)
nx= diff(xrng)/dx
ny= diff(yrng)/dy
targets= data.frame(
  x=rep(seq(from=min(xrng)+dx/2, to=max(xrng)-dx/2, length.out=nx), ny),
  y=sort(rep(seq(from=min(yrng)+dy/2, to=max(yrng)-dy/2, length.out=ny), nx))
)
plot(xrng,yrng,type="n",xlab="",ylab="", xaxt="n", yaxt="n" ,asp=1, bty="n")
clr= terrain.colors(ncolors)
# Show interpolation results (target locations)
for (i in 1:nrow(targets)) {
  weights= idweights(targets$x[i], targets$y[i], sources$x, sources$y,
    nsectors, norigins, power)
  if (abs(sum(weights)-1) > 0.0001) stop("Sum of weights not 1.")
  val= sum(weights * sources$value)
  clrind= round((val-min(sources$value))/
  (max(sources$value)-min(sources$value))*(ncolors-1))+1
  rect(xleft=targets$x[i]-dx,xright=targets$x[i]+dx,
    ybottom=targets$y[i]-dx,ytop=targets$y[i]+dx,col=clr[clrind],border=NA)
}
# Show original values (source locations)
vect_clrinds= round((sources$value-min(sources$value))/
  (max(sources$value)-min(sources$value))*(ncolors-1))+1
points(sources$x,sources$y,pch=20,cex=2,col=clr[vect_clrinds])
points(sources$x,sources$y,pch=1,cex=2)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

