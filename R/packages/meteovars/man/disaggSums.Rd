\name{disaggSums}
\alias{disaggSums}
\title{Disaggregation of sums based on a higher-resolution time series}
\usage{
  disaggSums(L, H, setZero = FALSE)
}
\arguments{
  \item{L}{A numeric vector representing the values of the
  lower-resolution series whose \emph{sum} is to be
  retained. The sum is only guaranteed to be retained
  exactly if argument \code{setZero} is \code{FALSE}.}

  \item{H}{A numeric vector representing the values of the
  higher-resolution series whose \emph{pattern} is to be
  retained. The lenght of \code{H} must be a multiple of
  the length of \code{L}, i.e. \code{length(H) %%
  length(L)} must evaluate to zero.}

  \item{setZero}{Logical value to control what happens if,
  for an element in \code{L}, the sum of the corresponding
  elements in \code{H} is zero.  If \code{setZero} is
  \code{FALSE} (default), the value of the element in
  \code{L} is 'equally' distributed over the elements of
  the result vector.  If \code{setZero} is \code{TRUE},
  however, those elements of the result vector are set to
  zero. See notes and examples.}
}
\value{
  A vector of the same lenght as \code{H} containing the
  result of the disaggregation.
}
\description{
  The function takes two vectors as input which usually
  represent the values of \emph{regular} time series. One
  of the vectors, \code{L}, is of low resolution and the
  other one, \code{L} is of higher resolution. Thus,
  \code{length(L)} is less than \code{length(H)}. This
  function creates a result vector where the \emph{sum} of
  the values is taken from \code{L} but the \emph{pattern}
  in terms of the variability in time is taken from
  \code{H}.
}
\note{
  A typical application is the disaggregation of rainfall.
  Using this method, it is possible, for example, to
  transfer the hourly pattern observed at an automatic rain
  gage to another observation with daily resolution made at
  a nearby rain gage. In that case, \code{setZero=FALSE} is
  a reasonable default, because a non-zero observation at
  the daily rain gage should be preserved, even if the
  nearby rain gage with hourly data did not record any
  precipitation in any 1-hour interval on that particular
  day.
}
\examples{
H= c(1,2,3,3,2,1,0,0,0,4,0,4)
L= c(2,4,6,8)
print(round(disaggSums(L, H, setZero=FALSE),3))
print(round(disaggSums(L, H, setZero=TRUE),3))
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

