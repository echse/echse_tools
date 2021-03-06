% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/topocatch.r
\name{dem.fill}
\alias{dem.fill}
\title{Filling of sinks in a digital elevation model (DEM)}
\usage{
dem.fill(fileIn, fileOut, ndigits, replace = FALSE, silent = TRUE)
}
\arguments{
\item{fileIn}{Name/path of input file (ASCII grid).}

\item{fileOut}{Name/path of output file (ASCII grid).}

\item{ndigits}{Number of relevant digits in the raw DEM (integer); see notes}

\item{replace}{Is is OK to replace an existing output file? Defaults to
  \code{FALSE}.}

\item{silent}{Print status info? (logical)}
}
\value{
\code{NULL}
}
\description{
See the workhorse function \code{\link{sinkfill}} for details.
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}

