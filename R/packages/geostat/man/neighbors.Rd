\name{neighbors}
\alias{neighbors}
\title{Find suitable neighboring points for spatial interpolation by sector search}
\usage{
  neighbors(x_tar, y_tar, x_src, y_src, nsectors, norigins)
}
\arguments{
  \item{x_tar}{x-coordinate of the target point (scalar)}

  \item{y_tar}{y-coordinate of the target point (scalar)}

  \item{x_src}{x-coordinates of the source points (vector)}

  \item{y_src}{y-coordinates of the source points (vector)}

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

  \item{nsectors}{The number of sectors to use (minimum 1).
  Each sector will be 360 / \code{nsectors} degrees wide.
  Setting \code{nsectors} to 1 results in a
  nearest-neighbor approach (also known as Thiessen
  method).}
}
\value{
  A list with the following two components:
  \item{selected}{A logical vector or the same length as
  \code{x_src}. The \code{TRUE} elements identify the
  selected source points. The order of elements is the same
  as in \code{x_src} and \code{y_src}.} \item{distance}{A
  numeric vector or the same length as \code{x_src}
  containing the distances between the target point and all
  source points.  The order of elements is the same as in
  \code{x_src} and \code{y_src}.}
}
\description{
  Given a point 'T' with coordinates \code{x_tar, y_tar}
  and a set of points 'S' defined by the coordinate vectors
  \code{x_src} and \code{y_src}, the function identifies
  those points among the set 'S' which are suitable for
  estimating the an unknown value of a variable at point
  'T' from the values at the points 'S'. This is basically
  the initial step of every spatial interpolation
  procedure. The function uses a sector search approach,
  i.e. the surrounding of point 'T' is divided into a
  user-specified number of sectors. For each sector, a
  \emph{single} point in the set 'S' is identified. This is
  the point in that particular sector which is closest to
  point 'T'. A sector can be empty.
}
\examples{
ndata= 16
nfigcol= 4
xrng= c(0, 100)
yrng= c(0, 100)
nsectors= 8
norigins= 10
data= data.frame(
  x= min(xrng) + runif(ndata) * diff(xrng),
  y= min(yrng) + runif(ndata) * diff(yrng)
)
old.par= par(mar=c(1,1,1,1))
split.screen(c(ceiling(ndata/nfigcol),nfigcol))
for (i in 1:nrow(data)) {
  x_tar= data$x[i]
  y_tar= data$y[i]
  x_src= data$x[1:nrow(data) != i]
  y_src= data$y[1:nrow(data) != i]
  nb= neighbors(x_tar, y_tar, x_src, y_src, nsectors, norigins)
  inds= which(nb$selected)
  avgdist= mean(nb$dist[nb$selected])
  screen(i)
  plot(data$x,data$y,type="n",xlab="x",ylab="y", xaxt="n", yaxt="n" ,asp=1)
  points(x_tar,y_tar,pch=1)
  points(x_src,y_src,pch=3,col="grey")
  for (k in 1:length(inds)) {
    lines(c(x_tar,x_src[inds[k]]), c(y_tar,y_src[inds[k]]), col="blue")
  }
  mtext(side=3,paste(round(avgdist,3)," (",sum(nb$selected),")",sep=""))
}
close.screen(all=TRUE)
par(old.par)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

