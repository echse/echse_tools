\name{geogrid.zones.continuous}
\alias{geogrid.zones.continuous}
\title{Zonal statistics for a grid of continuous data}
\usage{
  geogrid.zones.continuous(grid_zones, grid_data,
    minCoverage = 1, prefix = "", addMedian = FALSE,
    addQuartiles = FALSE, addExtremes = FALSE, sdigits = 6,
    silent = TRUE)
}
\arguments{
  \item{grid_zones}{A geogrid object defining the zones.
  The values in this grid should be integers.}

  \item{grid_data}{A geogrid object holding data of a
  continuous variable (like elevation, for example). The
  extend and resolution of this grid must be the same as
  for \code{grid_zones}.}

  \item{minCoverage}{A numeric value specifying the minimum
  coverage of the valid grid cells in \code{grid_zones} by
  the valid values in \code{grid_data}. The value must be
  in range \eqn{0 < minCoverage <= 1}.  If the actual
  coverage is less than the value of \code{minCoverage}, a
  warning is generated. Note that the coverage is analyzed
  globally but not for the individual zones.}

  \item{prefix}{A character string used as a prefix when
  creating column names for the output table. Must be
  compatible with R's convention for names.  See the return
  value for details.}

  \item{addMedian}{Should the median be added to the
  statistics? (logical)}

  \item{addQuartiles}{Should quartiles be added to the
  statistics? (logical)}

  \item{addExtremes}{Should min and max be added to the
  statistics? (logical)}

  \item{sdigits}{The number of significant digits in the
  output table. Should be an integer >= 1.}

  \item{silent}{Print diagnostic messages? (logical)}
}
\value{
  A data frame with at least two columns 'id' and 'avg'.
  The 'id' columns holds the zones' IDs (unique values in
  \code{grid_zone}). The 'avg' column contains the
  arithmetic mean of the spatially corresponding values in
  \code{grid_data}. If \code{addMedian} is \code{TRUE}, a
  column 'med' is appended. If \code{addQuartiles} is
  \code{TRUE}, the additional columns 'q25' and 'q75' are
  present which hold the values of the lower and upper
  quartiles.  If \code{addExtremes} is \code{TRUE}, there
  are two additional columns named 'min' and 'max'.  If the
  value of \code{prefix} is a non-empty string, this string
  is used as a prefix to the above-mentioned column names.
}
\description{
  The function expects two input grids, one with data of a
  continuous variable and another one with integer codes
  defining zones (classes). For each unique zone, the
  (spatially) corresponding values of the continuous
  variable are analyzed and a statistics table is returned.
  The function can be used, for example, to compute an
  elevation statistics for river basins.
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}

