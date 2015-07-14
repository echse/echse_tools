\name{externalInputLocationsTable}
\alias{externalInputLocationsTable}
\title{Generates an external input locations table for an ECHSE-based model}
\usage{
  externalInputLocationsTable(file_targets, files_sources,
    idField_targets = "id", idField_sources = "id",
    colsep = "\t", nsectors = 4, norigins = 4, power = 2,
    file_result, ndigits = 3, overwrite = FALSE)
}
\arguments{
  \item{file_targets}{File holding a table of target
  locations. The table must have (at least) three columns
  with names 'id' (default), 'x', and 'y'. The field
  separator can be specified with the \code{colsep}
  argument. The 'x' and 'y' field should contain numeric
  coordinates while the id-field contains unique
  identifiers for all target locations. In
  (semi)-distributed hydrological modeling, this is
  typically a table of the coordinates of the sub-basins'
  centers of gravity.}

  \item{files_sources}{A \emph{named} vector of file names.
  Each of the files is expected to hold a table of source
  locations for a particular variable.  The names of the
  variables are inferred from the names of the vectors'
  elements. In the context of hydrological modeling, these
  files typically list the coordinates of rain gages and
  sensors of other meteorological variables. All files must
  be in the same format as the file supplied as
  \code{file_targets}. The id-fields contain the station
  names/numbers.}

  \item{idField_targets}{Name of the field in
  \code{file_targets} containing ID strings for the
  locations. Defaults to 'id'.}

  \item{idField_sources}{Name of the field in
  \code{file_sources} containing ID strings for the
  locations. Defaults to 'id'.}

  \item{colsep}{Field separator used in all input files and
  the output file.  Defaults to TAB.}

  \item{nsectors}{See corresponding argument of the
  \code{\link{idweights}} method. The default value of 4 is
  often used and results in a so-called 'quadrant search'.}

  \item{norigins}{See corresponding argument of the
  \code{\link{idweights}} method.}

  \item{power}{See corresponding argument of the
  \code{\link{idweights}} method. Defaults to 2.}

  \item{file_result}{Name of the result file created by
  this function.}

  \item{ndigits}{The number of digits to be used when
  printing the computed weights to \code{file_result}. This
  value is passed as the second argument to the
  \code{\link{round}} method. Records for those source
  locations whose \emph{rounded} weight is zero are dropped
  when creating \code{file_result}.}

  \item{overwrite}{Is it OK to overwrite an existing result
  file? Defaults to FALSE.}
}
\value{
  \code{NULL}
}
\description{
  This method creates the so-called external input
  locations table for an ECHSE-based model. This table is
  used to assign a (set of) station(s) and weight(s) to all
  simulated objects for a particular variable. In
  hydrological models, the table usually assigns at least a
  (set of) raingage(s) to all sub-basins of a catchment.
}
\note{
  For more information on format and use of the result file
  please have a look at the ECHSE documentation (reference
  given below). In particular, one should read the
  information on external input variables and the
  configuration file's keyword
  'table_externalInput_locations'.
}
\examples{
# Demonstrates the assignment of hydro-meteorological sensors to sub-basins
# Dummy table of sub-basin coordinates
file_basins= tempfile()
tbl= data.frame(id= c("sb1","sb2"), x=runif(2), y=runif(2))
write.table(x= tbl, file=file_basins, sep="\\t", row.names=FALSE, quote=FALSE)
# Dummy table of rain gages
file_gages= tempfile()
tbl= data.frame(id= c("g1","g2","g3"), x=runif(3), y=runif(3))
write.table(x= tbl, file=file_gages, sep="\\t", row.names=FALSE, quote=FALSE)
# Dummy table of temperature sensors
file_temps= tempfile()
tbl= data.frame(id= c("t1","t2","t3","t4"), x=runif(4), y=runif(4))
write.table(x= tbl, file=file_temps, sep="\\t", row.names=FALSE, quote=FALSE)
# Compute and show result table
file_result= tempfile()
externalInputLocationsTable(
  file_targets= file_basins,
  files_sources= c(precip= file_gages, temper= file_temps),
  colsep="\\t", nsectors= 4, norigins= 2, power= 2,
  file_result=file_result)
print(read.table(file=file_result))
file.remove(file_basins, file_gages, file_temps, file_result)
}
\author{
  David Kneis \email{david.kneis@uni-potsdam.de}
}
\references{
  Kneis, D.: Eco-Hydrological Simulation Environment
  (ECHSE) - Documentation of the Generic Components.
}

