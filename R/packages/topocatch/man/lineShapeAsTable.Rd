% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/private_common_shapefile.r
\name{lineShapeAsTable}
\alias{lineShapeAsTable}
\title{Tabular representation of a shape file with line features}
\usage{
lineShapeAsTable(file, id_field, attribs = c(), endsOnly = FALSE)
}
\arguments{
\item{file}{Name/path of the shape file (with or without extension .shp). The
file must contain \emph{line} features (not points or polygons). The lines
must \emph{not} be multi-part lines, i.e. each line must be a continuous
sequence of coordinates without breaks (gaps).}

\item{id_field}{Name of the column in the shape file's attribute table
containing the unique object IDs.}

\item{attribs}{Vector of (further) attribute columns names. By default, no
further attributes are extracted.}

\item{endsOnly}{If \code{TRUE}, the result table contains only the
  coordinates of the lines' end points but no intermediate coordinates. If
  \code{FALSE} (default), the full array of coordinates is returned.}
}
\value{
A list of two data frames with names 'shp' and 'att'. The latter
  data frame represents the shape file's attribute table. It contains the
  ID field defined by \code{id_field} and the attribute fields listed in the
  \code{attribs} argument or \code{NULL} if \code{attribs} is an empty vector.
  The contents of the 'shp' data frame depends on the value of argument
  \code{endsOnly}. See notes for details.
}
\description{
The function converts an ESRI shape file with line features into a table of
coordinates. Feature attributes can be extracted as well.
}
\note{
If \code{endsOnly} is \code{FALSE}, the 'shp' data frame has the 3
  columns. The first column is the ID field with the name determined by
  \code{id_field}. The other two columns are 'x', and 'y'. The table
  contains the complete list of the lines' coordinates.
  The number of records in the result table is
  \eqn{\sum_{i=1}^n{L_i}} where \eqn{n} is the number of lines and \eqn{L_i}
  is the number of coordinates of a particular line with index \eqn{i}.

If \code{endsOnly} is \code{FALSE}, the 'shp' data frame has 5 columns.
The first column is the ID field with the name determined by \code{id_field}.
The 4 remaining columns are 'x1', 'y1', 'x2', and 'y2'. The table contains
only the coordinates of the lines' start and end points in the fields 'x1', 'y1' and
'x2', 'y2', respectively. The result table has as many records as there are
line features in the shape file.
}
\author{
David Kneis \email{david.kneis@uni-potsdam.de}
}
\seealso{
\code{\link{slot}}, \code{\link{readShapeLines}} from package \code{maptools}
}

