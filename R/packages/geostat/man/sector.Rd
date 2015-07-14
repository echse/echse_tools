\name{sector}
\alias{sector}
\title{Sector where a point is located when looking from another point}
\usage{
  sector(nsectors, angle_of_origin, x_from, y_from, x_to,
    y_to)
}
\arguments{
  \item{nsectors}{The number of sectors to use (minimum 1).
  Each sector will be 360 / \code{nsectors} degrees wide.
  Setting \code{nsectors} to 1 results in a
  nearest-neighbor approach (also known as Thiessen
  method).}

  \item{angle_of_origin}{Defines the angle of origin, i.e.
  the boundary between the first and the last sector. If
  \code{angle_of_origin} is zero, the boundary between
  sector 1 and sector \code{nsectors} is at zero degrees
  (i.e. the North-South axis). If, for example,
  \code{nsectors} is 4 and \code{angle_of_origin} is 90,
  the sector with index 1 covers the range from 90 to 180
  degrees.}

  \item{x_from}{x-coordinate of the base point}

  \item{y_from}{y-coordinate of the base point}

  \item{x_to}{x-coordinate of the point located from the
  base point}

  \item{y_to}{y-coordinate of the point located from the
  base point}
}
\value{
  The index of the sector as an integer in range 1 to
  \code{nsectors}.
}
\description{
  When looking from a point 'A' to another point 'B', the
  function computes the sector in which point 'B' is
  located. A sector is simply a range of angles. Sector
  indices increase in clock-wise direction. For example, if
  the number of sectors in 4, and the origin is at 0
  degrees, sector 1 would cover the range from 0 to 90
  degrees (North-East sector).
}
\examples{
nsources= 1000
xrng= c(0, 100)
yrng= c(0, 100)
x_tar= 50
y_tar= 50
nsectors= 4
orig_angle= 45.
src= data.frame(x= min(xrng) + runif(nsources) * diff(xrng),
  y= min(yrng) + runif(nsources) * diff(yrng))
src$sect= NA
for (i in 1:nrow(src))
  src$sect[i]= sector(nsectors,orig_angle,x_tar,y_tar,src$x[i],src$y[i])
plot(xrng, yrng, type="n", xlab="x", ylab="y", asp=1)
clr= rainbow(nsectors)
for (i in 1:nsectors)
  points(src$x[src$sect==i],src$y[src$sect==i],pch=i,cex=0.5,col=clr[i])
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

