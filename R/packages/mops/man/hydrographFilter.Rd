\name{hydrographFilter}
\alias{hydrographFilter}
\title{Filter a hydrograph for periods of high or low flow}
\usage{
hydrographFilter(vect, windowLength = ceiling(length(vect)/3),
  windowLag = ceiling(windowLength/10), alpha = 0.9, marginsOK = FALSE)
}
\arguments{
  \item{vect}{Numeric vector representing the hydrograph,
  possibly containing \code{NA}s (see notes). It is assumed
  that the values represent observations taken at equal
  time intervals (regular time series).}

  \item{windowLength}{Integer constant defining the width
  of the moving window (as number of observations) used in
  the data selection algorithm. Must be >= 1 and <= the
  length of \code{vect}. The window should be wide enough,
  i.e.  contain a reasonable number of values, to estimate
  the desired quantile.}

  \item{windowLag}{Integer constant defining the lag (as
  number of observations) being used when moving window
  along the time axis. Must be >= 1 and <= the length of
  \code{vect}. Typically, \code{windowLag} is chosen to be
  smaller than \code{windowLength} which results in an
  overlap of subsequent sampling windows. The computation
  time increases as \code{windowLag} decreases.}

  \item{alpha}{A numeric value between 0 and 1. A value of
  0.75 results in the selection of values exceeding the
  'within-window' 75 percentile, for example.}

  \item{marginsOK}{If \code{FALSE} (default), the algorithm
  only selects values in the center part of the moving
  window but not at the margins.}
}
\value{
A logical vector of the same length as \code{vect}.
Elements are \code{TRUE} where the corresponding values in
\code{vect} exceed the time- variable threshold defined by
\code{alpha}. See the examples for how this vector can be
used to filter \code{vect} for periods of floods or
droughts.
}
\description{
This method allows a hydrograph to be filtered for flood or
drought events. Flood (drought) events are those where the
observed variable exceeds (falls below) a threshold value.
The threshold value is conveniently defined as a quantile.
Furthermore, the threshold is variable in time. This is
achieved by applying the threshold-based selection concept
to a moving window being shifted along the time axis. The
quantile (alpha value) and the parameters of the moving
windows are user-specified.
}
\note{
Data will be selected only from those time windows where
all corresponding elements in \code{vect} are finite
according to \code{is.finite}. The selection algorithm
tests for exceedence of threshold using the > operator (not
>=). Thus, \code{TRUE} elements in the returned vector are
\emph{greater than} the threshold while \code{FALSE}
elements are \emph{equal or less than} the threshold.
}
\examples{
hydrograph= sin(seq(from=-pi/2, to=8*pi, by=0.05))+1
select_1= hydrographFilter(hydrograph, windowLag=1, alpha=0.8)
select_2= hydrographFilter(hydrograph, windowLag=1, alpha=0.2)
times= 1:length(hydrograph)
plot(x=times, y=hydrograph, xlab="Index", ylab="Data")
points(x=times[select_1],y=hydrograph[select_1], col="blue", pch=20)
points(x=times[!select_2],y=hydrograph[!select_2], col="red", pch=20)
legend("topright",pch=c(1,20,20),col=c("black","blue","red"),
  legend=c("Data","Floods","Droughts"))
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

